#!/usr/bin/env python3
"""Re-score saved grid models with a fine, restricted time-shift range."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from array import array
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from waveform_inversion import (
    COMPONENTS,
    MODEL_FILE,
    RESULTS,
    WORK,
    ensure_fk,
    parse_band,
    process_trace,
    rotate_rt_to_ne,
    station_records,
    synthesize,
    write_fk_model,
)


WINDOWS = (
    {"name": "pnl", "start": 0.0, "end": 90.0, "weight": 1.0},
    {"name": "surface", "start": 70.0, "end": 200.0, "weight": 1.0},
)


def read_grid(path: Path) -> list[dict]:
    with path.open() as f:
        return [
            {
                "depth": int(row["depth_km"]),
                "strike": int(row["strike"]),
                "dip": int(row["dip"]),
                "rake": int(row["rake"]),
                "old_cc": float(row["mean_windowed_cc"]),
            }
            for row in csv.DictReader(f)
        ]


def corr_at_shift(obs: array, syn: array, shift_samples: int) -> float | None:
    n = min(len(obs), len(syn))
    if shift_samples >= 0:
        o0, s0, m = shift_samples, 0, n - shift_samples
    else:
        o0, s0, m = 0, -shift_samples, n + shift_samples
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


def best_corr(obs: array, syn: array, shift_min: float, shift_max: float, shift_dt: float, sample_dt: float) -> tuple[float, float]:
    best = -2.0
    best_shift = shift_min
    n = int(round((shift_max - shift_min) / shift_dt))
    for i in range(n + 1):
        shift_s = shift_min + i * shift_dt
        shift_samples = int(round(shift_s / sample_dt))
        cc = corr_at_shift(obs, syn, shift_samples)
        if cc is not None and cc > best:
            best = cc
            best_shift = shift_samples * sample_dt
    return best, best_shift


def best_shared_corr(
    obs_by_comp: dict[str, array],
    syn_by_comp: dict[str, array],
    shift_min: float,
    shift_max: float,
    shift_dt: float,
    sample_dt: float,
) -> tuple[float, float, dict[str, float]]:
    best = -2.0
    best_shift = shift_min
    best_comp_ccs: dict[str, float] = {}
    n = int(round((shift_max - shift_min) / shift_dt))
    for i in range(n + 1):
        shift_s = shift_min + i * shift_dt
        shift_samples = int(round(shift_s / sample_dt))
        comp_ccs = {}
        for comp in COMPONENTS:
            cc = corr_at_shift(obs_by_comp[comp], syn_by_comp[comp], shift_samples)
            if cc is None:
                break
            comp_ccs[comp] = cc
        if len(comp_ccs) != len(COMPONENTS):
            continue
        mean_cc = sum(comp_ccs.values()) / len(COMPONENTS)
        if mean_cc > best:
            best = mean_cc
            best_shift = shift_samples * sample_dt
            best_comp_ccs = comp_ccs
    return best, best_shift, best_comp_ccs


def synthetic_by_component(sta: str, info: dict, model: dict, window: dict, dt: float, band: tuple[float, float], mw: float, duration: float) -> dict[str, array]:
    syn = synthesize(
        sta.replace(".", "_"),
        round(info["dist"], 1),
        info["az"],
        model["depth"],
        model["strike"],
        model["dip"],
        model["rake"],
        mw,
        duration,
    )
    sz = process_trace(syn["z"], window["start"], window["end"], dt, band)
    sr = process_trace(syn["r"], window["start"], window["end"], dt, band)
    st = process_trace(syn["t"], window["start"], window["end"], dt, band)
    sn, se = rotate_rt_to_ne(sr, st, info["az"])
    return {"BHZ": sz, "BH1": sn, "BH2": se}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(RESULTS / "grid_search_results.csv"))
    parser.add_argument("--sample-dt", type=float, default=0.1)
    parser.add_argument("--shift-min", type=float, default=2.0)
    parser.add_argument("--shift-max", type=float, default=6.0)
    parser.add_argument("--shift-dt", type=float, default=0.1)
    parser.add_argument("--bandpass", default="0.02/0.1")
    parser.add_argument("--mw", type=float, default=4.79)
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--limit", type=int, default=0, help="debug: only score the first N input models")
    args = parser.parse_args()

    band = parse_band(args.bandpass)
    if band is None:
        raise SystemExit("bandpass must not be off")

    models = read_grid(Path(args.input))
    if args.limit:
        models = models[: args.limit]
    stations = station_records()
    distances = [round(info["dist"], 1) for info in stations.values()]
    write_fk_model(MODEL_FILE, WORK / "Idaho_model")

    obs_cache: dict[tuple[str, str, str], array] = {}
    for sta, info in stations.items():
        for window in WINDOWS:
            for comp in COMPONENTS:
                obs_cache[(sta, window["name"], comp)] = process_trace(
                    info["components"][comp], window["start"], window["end"], args.sample_dt, band
                )

    rows = []
    details = []
    for imodel, model in enumerate(models, 1):
        ensure_fk(model["depth"], distances, 2048, 0.1)
        ccs = []
        shifts = []
        for sta, info in stations.items():
            for window in WINDOWS:
                syn = synthetic_by_component(sta, info, model, window, args.sample_dt, band, args.mw, args.duration)
                obs_by_comp = {comp: obs_cache[(sta, window["name"], comp)] for comp in COMPONENTS}
                mean_cc, shift, comp_ccs = best_shared_corr(
                    obs_by_comp,
                    syn,
                    args.shift_min,
                    args.shift_max,
                    args.shift_dt,
                    args.sample_dt,
                )
                ccs.append(mean_cc)
                shifts.append(shift)
                for comp, cc in comp_ccs.items():
                    details.append(
                        [
                            model["depth"],
                            model["strike"],
                            model["dip"],
                            model["rake"],
                            sta,
                            window["name"],
                            comp,
                            cc,
                            shift,
                        ]
                    )
        score = sum(ccs) / len(ccs)
        median_shift = sorted(shifts)[len(shifts) // 2]
        rows.append(
            [
                model["depth"],
                model["strike"],
                model["dip"],
                model["rake"],
                score,
                len(ccs),
                median_shift,
                model["old_cc"],
            ]
        )
        if imodel == 1 or imodel % 50 == 0:
            print(f"scored {imodel}/{len(models)} best_so_far={max(r[4] for r in rows):.4f}", flush=True)

    rows.sort(key=lambda row: row[4], reverse=True)
    tag = f"shift{args.shift_min:g}_{args.shift_max:g}_dt{args.shift_dt:g}".replace(".", "p")
    out_grid = RESULTS / f"grid_search_results_{tag}.csv"
    out_details = RESULTS / f"trace_scores_{tag}.csv"
    with out_grid.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "depth_km",
                "strike",
                "dip",
                "rake",
                "mean_windowed_cc",
                "n_window_traces",
                "median_shift_s",
                "previous_mean_windowed_cc",
            ]
        )
        writer.writerows(rows)
    with out_details.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["depth_km", "strike", "dip", "rake", "station", "window", "component", "cc", "shift_s"])
        writer.writerows(details)
    print(f"Wrote {out_grid}")
    print(f"Wrote {out_details}")


if __name__ == "__main__":
    main()
