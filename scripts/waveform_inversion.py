#!/usr/bin/env python3
"""Grid-search focal mechanism and depth with FK synthetics.

The script intentionally uses only the Python standard library so it can run in
the course VM without extra packages.  It reads SAC files directly, calls the
provided FK/syn programs, applies the long-period waveform-comparison steps
described in class, and scores each synthetic with time-shifted normalized
correlation.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import struct
import subprocess
from array import array
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FK_DIR = next(path for path in (ROOT / "sourceCode" / "fk", ROOT / "fk") if (path / "fk.pl").exists())
DISP_DIR = next(path for path in (ROOT / "data" / "disp", ROOT / "disp") if path.exists())
MODEL_FILE = next(path for path in (ROOT / "models" / "Idaho_Crustal_Model.txt", ROOT / "Idaho_Crustal_Model.txt") if path.exists())
RESULTS = ROOT / "results"
WORK = RESULTS / "work"

COMPONENTS = ("BHZ", "BH1", "BH2")
DEFAULT_WINDOWS = "pnl:0:90:1.0,surface:70:200:1.0"


def read_sac(path: Path) -> dict:
    raw = path.read_bytes()
    floats = struct.unpack("<70f", raw[:280])
    ints = struct.unpack("<40i", raw[280:440])
    chars = raw[440:632]
    npts = ints[9]
    data = array("f")
    data.frombytes(raw[632 : 632 + 4 * npts])
    if data.itemsize != 4:
        data.byteswap()
    keys = [
        chars[i * 8 : (i + 1) * 8].decode("ascii", "ignore").strip("\x00 ")
        for i in range(24)
    ]
    return {
        "path": path,
        "data": data,
        "delta": floats[0],
        "b": floats[5],
        "e": floats[6],
        "dist": floats[50],
        "az": floats[51],
        "baz": floats[52],
        "cmpaz": floats[57],
        "cmpinc": floats[58],
        "npts": npts,
        "kstnm": keys[0],
        "kcmpnm": keys[20],
    }


def write_fk_model(source: Path, target: Path) -> None:
    lines = []
    for line in source.read_text().splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        vals = text.split()
        if len(vals) < 6:
            continue
        vp, vs, rho, thick, qp, qs = vals[:6]
        lines.append(f"{float(thick):.4f} {float(vs):.4f} {float(vp):.4f} {float(rho):.4f} {float(qs):.1f} {float(qp):.1f}\n")
    target.write_text("".join(lines))


def station_records() -> dict:
    stations: dict[str, dict] = {}
    for path in sorted(DISP_DIR.glob("*.sac.dis")):
        rec = read_sac(path)
        sta = f"{path.name.split('.')[0]}.{path.name.split('.')[1]}"
        stations.setdefault(sta, {"components": {}})
        stations[sta]["components"][rec["kcmpnm"]] = rec
        stations[sta]["dist"] = rec["dist"]
        stations[sta]["az"] = rec["az"]
    return stations


def resample_to_grid(rec: dict, start: float, end: float, dt: float) -> array:
    n = int(round((end - start) / dt)) + 1
    out = array("f", [0.0]) * n
    data = rec["data"]
    b = rec["b"]
    odt = rec["delta"]
    last = rec["npts"] - 1
    for i in range(n):
        t = start + i * dt
        x = (t - b) / odt
        j = int(math.floor(x))
        if j < 0 or j >= last:
            continue
        frac = x - j
        out[i] = data[j] * (1.0 - frac) + data[j + 1] * frac
    demean(out)
    return out


def demean(values: array) -> None:
    if not values:
        return
    mean = sum(values) / len(values)
    for i, val in enumerate(values):
        values[i] = val - mean


def biquad_filter(values: array, dt: float, cutoff: float, kind: str) -> array:
    """Apply one RBJ cookbook low-pass or high-pass biquad."""
    if cutoff <= 0.0:
        return array("f", values)
    nyquist = 0.5 / dt
    cutoff = min(cutoff, nyquist * 0.98)
    w0 = 2.0 * math.pi * cutoff * dt
    cosw = math.cos(w0)
    sinw = math.sin(w0)
    q = math.sqrt(0.5)
    alpha = sinw / (2.0 * q)
    if kind == "lowpass":
        b0 = (1.0 - cosw) / 2.0
        b1 = 1.0 - cosw
        b2 = (1.0 - cosw) / 2.0
    elif kind == "highpass":
        b0 = (1.0 + cosw) / 2.0
        b1 = -(1.0 + cosw)
        b2 = (1.0 + cosw) / 2.0
    else:
        raise ValueError(f"unknown filter kind: {kind}")
    a0 = 1.0 + alpha
    a1 = -2.0 * cosw
    a2 = 1.0 - alpha
    b0, b1, b2, a1, a2 = b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0
    out = array("f", [0.0]) * len(values)
    x1 = x2 = y1 = y2 = 0.0
    for i, x0 in enumerate(values):
        y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2
        out[i] = y0
        x2, x1 = x1, x0
        y2, y1 = y1, y0
    return out


def zero_phase_bandpass(values: array, dt: float, fmin: float, fmax: float, passes: int = 2) -> array:
    """Small dependency-free long-period band-pass used before correlation."""
    out = array("f", values)
    for _ in range(max(1, passes)):
        out = biquad_filter(out, dt, fmin, "highpass")
        out = biquad_filter(out, dt, fmax, "lowpass")
    out.reverse()
    for _ in range(max(1, passes)):
        out = biquad_filter(out, dt, fmin, "highpass")
        out = biquad_filter(out, dt, fmax, "lowpass")
    out.reverse()
    demean(out)
    return out


def process_trace(rec: dict, start: float, end: float, dt: float, band: tuple[float, float] | None) -> array:
    trace = resample_to_grid(rec, start, end, dt)
    if band is not None:
        trace = zero_phase_bandpass(trace, dt, band[0], band[1])
    return trace


def corr_with_shift(obs: array, syn: array, max_shift: int) -> tuple[float, int]:
    best = -2.0
    best_shift = 0
    n = min(len(obs), len(syn))
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            o0, s0, m = shift, 0, n - shift
        else:
            o0, s0, m = 0, -shift, n + shift
        if m < 20:
            continue
        dot = oo = ss = 0.0
        for k in range(m):
            o = obs[o0 + k]
            s = syn[s0 + k]
            dot += o * s
            oo += o * o
            ss += s * s
        if oo <= 0.0 or ss <= 0.0:
            continue
        cc = dot / math.sqrt(oo * ss)
        if cc > best:
            best = cc
            best_shift = shift
    return best, best_shift


def corr_at_shift(obs: array, syn: array, shift: int) -> float | None:
    n = min(len(obs), len(syn))
    if shift >= 0:
        o0, s0, m = shift, 0, n - shift
    else:
        o0, s0, m = 0, -shift, n + shift
    if m < 20:
        return None
    dot = oo = ss = 0.0
    for k in range(m):
        o = obs[o0 + k]
        s = syn[s0 + k]
        dot += o * s
        oo += o * o
        ss += s * s
    if oo <= 0.0 or ss <= 0.0:
        return None
    return dot / math.sqrt(oo * ss)


def corr_shared_shift(obs_by_comp: dict[str, array], syn_by_comp: dict[str, array], max_shift: int) -> tuple[float, int, dict[str, float]]:
    """Find one shift shared by all three components in one station/window."""
    best = -2.0
    best_shift = 0
    best_comp_ccs: dict[str, float] = {}
    for shift in range(-max_shift, max_shift + 1):
        comp_ccs = {}
        for comp in COMPONENTS:
            cc = corr_at_shift(obs_by_comp[comp], syn_by_comp[comp], shift)
            if cc is None:
                break
            comp_ccs[comp] = cc
        if len(comp_ccs) != len(COMPONENTS):
            continue
        mean_cc = sum(comp_ccs.values()) / len(COMPONENTS)
        if mean_cc > best:
            best = mean_cc
            best_shift = shift
            best_comp_ccs = comp_ccs
    return best, best_shift, best_comp_ccs


def parse_windows(text: str) -> list[dict]:
    windows = []
    for item in text.split(","):
        parts = item.split(":")
        if len(parts) not in (3, 4):
            raise ValueError(f"bad window specification: {item}")
        name, start, end = parts[:3]
        weight = float(parts[3]) if len(parts) == 4 else 1.0
        windows.append({"name": name, "start": float(start), "end": float(end), "weight": weight})
    if not windows:
        raise ValueError("at least one scoring window is required")
    return windows


def parse_band(text: str) -> tuple[float, float] | None:
    if text.lower() in {"none", "off", "0"}:
        return None
    fmin, fmax = (float(x) for x in text.split("/", 1))
    if not 0.0 < fmin < fmax:
        raise ValueError("bandpass must satisfy 0 < fmin < fmax")
    return fmin, fmax


def run(cmd: list[str], cwd: Path) -> None:
    env = os.environ.copy()
    env["PATH"] = f"{FK_DIR}:{env.get('PATH', '')}"
    subprocess.run(cmd, cwd=cwd, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def ensure_fk(depth: int, distances: list[float], nt: int, dt: float, force: bool = False) -> Path:
    model = WORK / "Idaho_model"
    out_dir = WORK / f"Idaho_model_{depth:g}"
    needed = [out_dir / f"{dist:.1f}.fk.0" for dist in distances]
    if force and out_dir.exists():
        shutil.rmtree(out_dir)
    if all(path.exists() for path in needed):
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "perl",
        str(FK_DIR / "fk.pl"),
        f"-M{model.name}/{depth}",
        f"-N{nt}/{dt}/1/0.2",
    ] + [f"{d:.1f}" for d in distances]
    run(cmd, WORK)
    return out_dir


def synthesize(sta: str, dist: float, az: float, depth: int, strike: int, dip: int, rake: int, mw: float, dura: float, force: bool = False) -> dict:
    out_base = WORK / "syn" / f"{sta}_{depth}_{strike}_{dip}_{rake}.z"
    out_base.parent.mkdir(parents=True, exist_ok=True)
    gf = WORK / f"Idaho_model_{depth:g}" / f"{dist:.1f}.fk.0"
    cmd = [
        str(FK_DIR / "syn"),
        f"-M{mw:.2f}/{strike}/{dip}/{rake}",
        "-I",
        f"-D{dura:g}",
        f"-A{az:.3f}",
        f"-O{out_base.name}",
        f"-G{gf}",
    ]
    if force or not (out_base.exists() and out_base.with_suffix(".r").exists() and out_base.with_suffix(".t").exists()):
        run(cmd, out_base.parent)
    z = read_sac(out_base)
    r = read_sac(out_base.with_suffix(".r"))
    t = read_sac(out_base.with_suffix(".t"))
    return {"z": z, "r": r, "t": t}


def rotate_rt_to_ne(r: array, t: array, az: float) -> tuple[array, array]:
    rad = math.radians(az)
    ca = math.cos(rad)
    sa = math.sin(rad)
    north = array("f", [0.0]) * len(r)
    east = array("f", [0.0]) * len(r)
    for i in range(len(r)):
        north[i] = r[i] * ca - t[i] * sa
        east[i] = r[i] * sa + t[i] * ca
    return north, east


def frange_int(start: int, stop: int, step: int) -> list[int]:
    return list(range(start, stop + (1 if step > 0 else -1), step))


def search(args: argparse.Namespace) -> None:
    RESULTS.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    write_fk_model(MODEL_FILE, WORK / "Idaho_model")
    stations = station_records()
    distances = [round(info["dist"], 1) for info in stations.values()]
    windows = parse_windows(args.windows)
    band = parse_band(args.bandpass)

    obs_grid: dict[str, dict[str, dict[str, array]]] = {}
    for sta, info in stations.items():
        obs_grid[sta] = {}
        for window in windows:
            obs_grid[sta][window["name"]] = {}
            for comp in COMPONENTS:
                obs_grid[sta][window["name"]][comp] = process_trace(
                    info["components"][comp], window["start"], window["end"], args.score_dt, band
                )

    rows = []
    details = []
    best = None
    max_shift = int(round(args.max_shift / args.score_dt))

    strikes = frange_int(args.strike_min, args.strike_max, args.strike_step)
    dips = frange_int(args.dip_min, args.dip_max, args.dip_step)
    rakes = frange_int(args.rake_min, args.rake_max, args.rake_step)
    depths = frange_int(args.depth_min, args.depth_max, args.depth_step)

    for depth in depths:
        ensure_fk(depth, distances, args.fk_nt, args.fk_dt, args.force_fk)
        for strike in strikes:
            for dip in dips:
                for rake in rakes:
                    ccs = []
                    shifts = []
                    for sta, info in stations.items():
                        syn = synthesize(
                            sta.replace(".", "_"), round(info["dist"], 1), info["az"], depth,
                            strike, dip, rake, args.mw, args.duration, args.force_syn
                        )
                        for window in windows:
                            sz = process_trace(syn["z"], window["start"], window["end"], args.score_dt, band)
                            sr = process_trace(syn["r"], window["start"], window["end"], args.score_dt, band)
                            st = process_trace(syn["t"], window["start"], window["end"], args.score_dt, band)
                            sn, se = rotate_rt_to_ne(sr, st, info["az"])
                            syn_by_comp = {"BHZ": sz, "BH1": sn, "BH2": se}
                            obs_by_comp = obs_grid[sta][window["name"]]
                            mean_cc, sh, comp_ccs = corr_shared_shift(obs_by_comp, syn_by_comp, max_shift)
                            if mean_cc > -1.5:
                                ccs.append(mean_cc * window["weight"])
                                shifts.append(sh * args.score_dt)
                                if args.keep_details:
                                    for comp, cc in comp_ccs.items():
                                        details.append([depth, strike, dip, rake, sta, window["name"], comp, cc, sh * args.score_dt])
                    score = sum(ccs) / sum(window["weight"] * len(stations) for window in windows)
                    med_shift = sorted(shifts)[len(shifts) // 2] if shifts else 0.0
                    row = [depth, strike, dip, rake, score, len(ccs), med_shift]
                    rows.append(row)
                    if best is None or score > best[4]:
                        best = row
                        print(
                            f"best cc={score:.4f} depth={depth} strike={strike} dip={dip} rake={rake} traces={len(ccs)}",
                            flush=True,
                        )

    rows.sort(key=lambda x: x[4], reverse=True)
    with (RESULTS / "grid_search_results.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["depth_km", "strike", "dip", "rake", "mean_windowed_cc", "n_window_traces", "median_shift_s"])
        writer.writerows(rows)
    if details:
        with (RESULTS / "trace_scores.csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["depth_km", "strike", "dip", "rake", "station", "window", "component", "cc", "shift_s"])
            writer.writerows(details)
    write_summary(rows, stations, args)


def write_summary(rows: list[list], stations: dict, args: argparse.Namespace) -> None:
    best = rows[0]
    lines = [
        "# Waveform Inversion Result",
        "",
        "## Data and Method",
        "",
        f"- Observed data: {len(stations)} stations x 3 components, instrument-response-removed displacement SAC files in `disp/`.",
        "- Synthetics: FK Green functions plus `syn`, using the supplied Idaho layered crustal model.",
        f"- Scoring windows: `{args.windows}` s relative to origin time; resampled to {args.score_dt:g} s.",
        f"- Long-period band-pass: `{args.bandpass}` Hz before waveform comparison.",
        f"- Misfit metric: maximum normalized correlation over +/-{args.max_shift:g} s time shifts in each window; objective is weighted mean correlation.",
        "",
        "## Best Model",
        "",
        f"- Depth: {best[0]} km",
        f"- Strike: {best[1]} deg",
        f"- Dip: {best[2]} deg",
        f"- Rake: {best[3]} deg",
        f"- Mean shifted windowed correlation: {best[4]:.4f} over {best[5]} windowed traces",
        "",
        "## Top 10 Grid Points",
        "",
        "| Rank | Depth (km) | Strike | Dip | Rake | Mean windowed CC | Median shift (s) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for i, row in enumerate(rows[:10], 1):
        lines.append(f"| {i} | {row[0]} | {row[1]} | {row[2]} | {row[3]} | {row[4]:.4f} | {row[6]:.1f} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "This run follows the class workflow more closely than the first pass: the traces are long-period filtered and scored in separated Pnl/surface-wave windows. "
        "Residual time shifts remain because the one-dimensional crustal model is still an approximation and no station corrections are solved for.",
    ]
    (RESULTS / "waveform_inversion_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth-min", type=int, default=6)
    parser.add_argument("--depth-max", type=int, default=18)
    parser.add_argument("--depth-step", type=int, default=4)
    parser.add_argument("--strike-min", type=int, default=0)
    parser.add_argument("--strike-max", type=int, default=330)
    parser.add_argument("--strike-step", type=int, default=30)
    parser.add_argument("--dip-min", type=int, default=30)
    parser.add_argument("--dip-max", type=int, default=90)
    parser.add_argument("--dip-step", type=int, default=30)
    parser.add_argument("--rake-min", type=int, default=-180)
    parser.add_argument("--rake-max", type=int, default=180)
    parser.add_argument("--rake-step", type=int, default=60)
    parser.add_argument("--mw", type=float, default=4.79)
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--windows", default=DEFAULT_WINDOWS, help="comma-separated name:start:end:weight windows")
    parser.add_argument("--bandpass", default="0.02/0.1", help="long-period bandpass in Hz, or off")
    parser.add_argument("--score-dt", type=float, default=0.5)
    parser.add_argument("--max-shift", type=float, default=40.0)
    parser.add_argument("--fk-nt", type=int, default=2048)
    parser.add_argument("--fk-dt", type=float, default=0.1)
    parser.add_argument("--force-fk", action="store_true")
    parser.add_argument("--force-syn", action="store_true")
    parser.add_argument("--keep-details", action="store_true")
    args = parser.parse_args()
    search(args)


if __name__ == "__main__":
    main()
