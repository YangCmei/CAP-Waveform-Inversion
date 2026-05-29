#!/usr/bin/env python3
"""Plot observed vs synthetic waveforms for the top grid-search models."""

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


COMP_LABELS = {
    "BHZ": "BHZ / Z",
    "BH1": "BH1 / N",
    "BH2": "BH2 / E",
}


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 12,
    anchor: str = "middle",
    color: str = "#20242a",
    weight: str = "400",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Microsoft YaHei, sans-serif" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{color}">{text}</text>'
    )


def top_models(limit: int, path: Path) -> list[dict]:
    with path.open() as f:
        rows = list(csv.DictReader(f))
    out = []
    for row in rows[:limit]:
        out.append(
            {
                "depth": int(row["depth_km"]),
                "strike": int(row["strike"]),
                "dip": int(row["dip"]),
                "rake": int(row["rake"]),
                "cc": float(row["mean_windowed_cc"]),
                "median_shift": float(row["median_shift_s"]),
            }
        )
    return out


def refined_summary() -> dict[tuple[int, int, int, int], dict]:
    path = RESULTS / "refined_top_model_scores.csv"
    if not path.exists():
        return {}
    out = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (int(row["depth_km"]), int(row["strike"]), int(row["dip"]), int(row["rake"]))
            out[key] = {
                "fine_cc": float(row["fine_mean_windowed_cc"]),
                "fine_median_shift": float(row["fine_median_shift_s"]),
            }
    return out


def refined_plot_shifts() -> dict[tuple[int, int, int, int, str, str], float]:
    path = RESULTS / "refined_trace_scores_top5.csv"
    if not path.exists():
        return {}
    grouped: dict[tuple[int, int, int, int, str, str], list[float]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (
                int(row["depth_km"]),
                int(row["strike"]),
                int(row["dip"]),
                int(row["rake"]),
                row["station"],
                row["component"],
            )
            grouped.setdefault(key, []).append(float(row["fine_shift_s"]))
    return {key: sorted(vals)[len(vals) // 2] for key, vals in grouped.items()}


def station_plot_shifts(path: Path) -> dict[tuple[int, int, int, int, str], float]:
    if not path.exists():
        return {}
    grouped: dict[tuple[int, int, int, int, str], list[float]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (
                int(row["depth_km"]),
                int(row["strike"]),
                int(row["dip"]),
                int(row["rake"]),
                row["station"],
            )
            grouped.setdefault(key, []).append(float(row["shift_s"]))
    return {key: sorted(vals)[len(vals) // 2] for key, vals in grouped.items()}


def normalize_trace(values: array) -> list[float]:
    scale = max((abs(v) for v in values), default=0.0)
    if scale <= 0.0 or not math.isfinite(scale):
        return [0.0 for _ in values]
    return [v / scale for v in values]


def path_for_trace(
    values: list[float],
    start: float,
    dt: float,
    plot_start: float,
    plot_end: float,
    x0: float,
    y0: float,
    width: float,
    amp: float,
    time_shift: float = 0.0,
    max_points: int = 500,
) -> str:
    if not values:
        return ""
    step = max(1, len(values) // max_points)
    commands = []
    drawing = False
    for i in range(0, len(values), step):
        t = start + i * dt + time_shift
        if t < plot_start or t > plot_end:
            drawing = False
            continue
        x = x0 + (t - plot_start) / (plot_end - plot_start) * width
        y = y0 - values[i] * amp
        commands.append(("L" if drawing else "M") + f" {x:.1f} {y:.1f}")
        drawing = True
    return " ".join(commands)


def station_order(stations: dict) -> list[str]:
    return sorted(stations, key=lambda name: stations[name]["dist"])


def synthetic_components(sta: str, info: dict, model: dict, start: float, end: float, dt: float, band: tuple[float, float], mw: float, duration: float) -> dict[str, array]:
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
    sz = process_trace(syn["z"], start, end, dt, band)
    sr = process_trace(syn["r"], start, end, dt, band)
    st = process_trace(syn["t"], start, end, dt, band)
    sn, se = rotate_rt_to_ne(sr, st, info["az"])
    return {"BHZ": sz, "BH1": sn, "BH2": se}


def make_svg(
    model: dict,
    rank: int,
    stations: dict,
    start: float,
    end: float,
    dt: float,
    band: tuple[float, float],
    mw: float,
    duration: float,
    apply_shift: bool,
    station_shifts: dict[tuple[int, int, int, int, str], float],
) -> str:
    ordered = station_order(stations)
    width, height = 1800, 1360
    left, right, top, bottom = 190, 45, 155, 90
    panel_gap = 38
    panel_w = (width - left - right - 2 * panel_gap) / 3
    row_h = (height - top - bottom) / len(ordered)
    amp = row_h * 0.34
    key = (model["depth"], model["strike"], model["dip"], model["rake"])
    score = model["cc"]
    shift = model["median_shift"] if apply_shift else 0.0

    title = (
        f"Rank {rank}: Depth {model['depth']} km, Strike {model['strike']}°, "
        f"Dip {model['dip']}°, Rake {model['rake']}°, CC={score:.4f}"
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 44, title, 28, weight="700"),
        f'<rect x="{width - 405}" y="74" width="345" height="46" rx="4" fill="#ffffff" stroke="#d8dde3"/>',
        f'<line x1="{width - 382}" y1="93" x2="{width - 332}" y2="93" stroke="#9aa3ad" stroke-width="2.4"/>',
        svg_text(width - 322, 98, "Observed", 15, "start", "#30343b", "700"),
        f'<line x1="{width - 190}" y1="93" x2="{width - 140}" y2="93" stroke="#1f6f8b" stroke-width="2.6"/>',
        svg_text(width - 130, 98, "Synthetic", 15, "start", "#1f6f8b", "700"),
    ]

    obs_grid: dict[str, dict[str, array]] = {}
    syn_grid: dict[str, dict[str, array]] = {}
    distances = [round(info["dist"], 1) for info in stations.values()]
    write_fk_model(MODEL_FILE, WORK / "Idaho_model")
    ensure_fk(model["depth"], distances, 2048, 0.1)
    for sta in ordered:
        info = stations[sta]
        obs_grid[sta] = {
            comp: process_trace(info["components"][comp], start, end, dt, band)
            for comp in COMPONENTS
        }
        syn_grid[sta] = synthetic_components(sta, info, model, start, end, dt, band, mw, duration)

    for icomp, comp in enumerate(COMPONENTS):
        x0 = left + icomp * (panel_w + panel_gap)
        parts += [
            f'<rect x="{x0:.1f}" y="{top - 22:.1f}" width="{panel_w:.1f}" height="{height - top - bottom + 28:.1f}" fill="#f7f8fa" stroke="#d8dde3"/>',
            svg_text(x0 + panel_w / 2, top - 36, COMP_LABELS[comp], 19, weight="700"),
            f'<line x1="{x0:.1f}" y1="{height - bottom + 8:.1f}" x2="{x0 + panel_w:.1f}" y2="{height - bottom + 8:.1f}" stroke="#6b7280" stroke-width="1.1"/>',
            svg_text(x0, height - bottom + 33, f"{start:g}s", 15, "start", "#20242a", "700"),
            svg_text(x0 + panel_w, height - bottom + 33, f"{end:g}s", 15, "end", "#20242a", "700"),
        ]
        for tick in range(int(start), int(end) + 1, 50):
            x = x0 + (tick - start) / (end - start) * panel_w
            parts.append(f'<line x1="{x:.1f}" y1="{top - 10:.1f}" x2="{x:.1f}" y2="{height - bottom + 8:.1f}" stroke="#e5e7eb" stroke-width="1"/>')

        for ista, sta in enumerate(ordered):
            y = top + ista * row_h + row_h / 2
            obs = normalize_trace(obs_grid[sta][comp])
            syn = normalize_trace(syn_grid[sta][comp])
            obs_path = path_for_trace(obs, start, dt, start, end, x0, y, panel_w, amp)
            trace_shift = station_shifts.get(key + (sta,), shift) if apply_shift else 0.0
            syn_path = path_for_trace(syn, start, dt, start, end, x0, y, panel_w, amp, time_shift=trace_shift)
            parts.append(f'<path d="{obs_path}" fill="none" stroke="#9aa3ad" stroke-width="0.9"/>')
            parts.append(f'<path d="{syn_path}" fill="none" stroke="#1f6f8b" stroke-width="1.25"/>')
            if icomp == 0:
                parts.append(svg_text(x0 - 12, y + 4, f'{sta} ({stations[sta]["dist"]:.0f} km)', 13, "end", "#20242a", "700"))

    parts += [
        svg_text(width / 2, height - 20, "Time (s)", 18, color="#20242a", weight="700"),
        "</svg>",
    ]
    return "\n".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--bandpass", default="0.02/0.1")
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--start", type=float, default=0.0)
    parser.add_argument("--end", type=float, default=200.0)
    parser.add_argument("--mw", type=float, default=4.79)
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--no-shift", action="store_true", help="plot synthetics without median-shift alignment")
    parser.add_argument(
        "--results-file",
        default=str(RESULTS / "grid_search_results_shift2_6_dt0p1.csv")
        if (RESULTS / "grid_search_results_shift2_6_dt0p1.csv").exists()
        else str(RESULTS / "grid_search_results.csv"),
    )
    parser.add_argument(
        "--trace-scores-file",
        default=str(RESULTS / "trace_scores_shift2_6_dt0p1.csv")
        if (RESULTS / "trace_scores_shift2_6_dt0p1.csv").exists()
        else "",
    )
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    stations = station_records()
    band = parse_band(args.bandpass)
    if band is None:
        raise SystemExit("bandpass must not be off for this plot")

    index_lines = [
        "# Top waveform-fit figures",
        "",
        "Gray = observed low-frequency waveform; blue = synthetic waveform.",
        "",
    ]
    station_shifts = station_plot_shifts(Path(args.trace_scores_file)) if args.trace_scores_file else {}
    for rank, model in enumerate(top_models(args.top, Path(args.results_file)), 1):
        svg = make_svg(
            model,
            rank,
            stations,
            args.start,
            args.end,
            args.dt,
            band,
            args.mw,
            args.duration,
            not args.no_shift,
            station_shifts,
        )
        name = (
            f"waveform_fit_rank{rank:02d}_D{model['depth']}_S{model['strike']}_"
            f"Dip{model['dip']}_R{model['rake']}.svg"
        )
        path = RESULTS / name
        path.write_text(svg)
        index_lines.append(
            f"- Rank {rank}: `{name}` "
            f"Depth={model['depth']} km, Strike={model['strike']}°, Dip={model['dip']}°, "
            f"Rake={model['rake']}°, CC={model['cc']:.4f}"
        )
        print(f"Wrote {path}")
    (RESULTS / "waveform_fit_top5_index.md").write_text("\n".join(index_lines) + "\n")
    print(f"Wrote {RESULTS / 'waveform_fit_top5_index.md'}")


if __name__ == "__main__":
    main()
