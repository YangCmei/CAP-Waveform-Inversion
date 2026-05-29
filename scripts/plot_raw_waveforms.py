#!/usr/bin/env python3
"""Plot raw instrument-response-removed displacement waveforms."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from waveform_inversion import COMPONENTS, RESULTS, station_records


COMP_LABELS = {
    "BH1": "BH1 / N",
    "BH2": "BH2 / E",
    "BHZ": "BHZ / Z",
}


def svg_text(x: float, y: float, text: str, size: int = 12, anchor: str = "middle", color: str = "#20242a") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Microsoft YaHei, sans-serif" text-anchor="{anchor}" '
        f'fill="{color}">{text}</text>'
    )


def normalize(values: list[float | None]) -> list[float | None]:
    scale = max((abs(v) for v in values if v is not None), default=0.0)
    if scale <= 0.0 or not math.isfinite(scale):
        return [None if v is None else 0.0 for v in values]
    return [None if v is None else v / scale for v in values]


def raw_trace_on_grid(rec: dict, start: float, end: float, dt: float) -> list[float | None]:
    n = int(round((end - start) / dt)) + 1
    out: list[float | None] = [None] * n
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
    return out


def trace_path(values: list[float | None], start: float, end: float, x0: float, y0: float, width: float, amp: float, max_points: int = 1200) -> str:
    if not values:
        return ""
    n = len(values)
    step = max(1, n // max_points)
    denom = max(1, n - 1)
    commands = []
    drawing = False
    for i in range(0, n, step):
        if values[i] is None:
            drawing = False
            continue
        t = start + (end - start) * i / denom
        x = x0 + (t - start) / (end - start) * width
        y = y0 - values[i] * amp  # type: ignore[operator]
        commands.append(("L" if drawing else "M") + f" {x:.1f} {y:.1f}")
        drawing = True
    return " ".join(commands)


def ordered_stations(stations: dict, sort_by: str) -> list[str]:
    if sort_by == "az":
        return sorted(stations, key=lambda name: stations[name]["az"])
    if sort_by == "name":
        return sorted(stations)
    return sorted(stations, key=lambda name: stations[name]["dist"])


def station_axis_label(stations: dict, sta: str, sort_by: str) -> str:
    if sort_by == "az":
        return f'{sta} ({stations[sta]["az"]:.0f} deg)'
    if sort_by == "dist":
        return f'{sta} ({stations[sta]["dist"]:.0f} km)'
    return sta


def default_time_range(stations: dict) -> tuple[float, float]:
    records = [
        rec
        for info in stations.values()
        for rec in info["components"].values()
    ]
    data_start = min(rec["b"] for rec in records)
    data_end = max(rec["e"] for rec in records)
    return math.floor(data_start / 200.0) * 200.0, math.ceil(data_end / 200.0) * 200.0


def time_suffix(start: float, end: float, is_default_window: bool) -> str:
    if is_default_window:
        return ""
    return f"_{start:g}_{end:g}s".replace(".", "p")


def tick_step(start: float, end: float) -> int:
    duration = end - start
    if duration <= 250:
        return 50
    if duration <= 800:
        return 100
    return 200


def time_ticks(start: float, end: float) -> list[int]:
    step = tick_step(start, end)
    first = math.ceil(start / step) * step
    return list(range(int(first), int(end) + 1, step))


def make_record_section(stations: dict, start: float, end: float, dt: float, sort_by: str) -> str:
    ordered = ordered_stations(stations, sort_by)
    width, height = 1500, 930
    margin_left, margin_right = 170, 40
    top, bottom = 135, 65
    panel_gap = 36
    panel_w = (width - margin_left - margin_right - 2 * panel_gap) / 3
    row_h = (height - top - bottom) / len(ordered)
    amp = row_h * 0.36

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 42, "Normalized displacement record sections", 28),
    ]

    for icomp, comp in enumerate(COMPONENTS):
        x0 = margin_left + icomp * (panel_w + panel_gap)
        parts += [
            f'<rect x="{x0:.1f}" y="{top - 24:.1f}" width="{panel_w:.1f}" height="{height - top - bottom + 28:.1f}" fill="#f7f8fa" stroke="#d8dde3"/>',
            svg_text(x0 + panel_w / 2, top - 38, COMP_LABELS.get(comp, comp), 18),
        ]
        for tick in time_ticks(start, end):
            x = x0 + (tick - start) / (end - start) * panel_w
            parts.append(f'<line x1="{x:.1f}" y1="{top - 12:.1f}" x2="{x:.1f}" y2="{height - bottom + 8:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
            parts.append(svg_text(x, height - bottom + 30, f"{tick:g}", 13, "middle", "#20242a"))
        for ista, sta in enumerate(ordered):
            rec = stations[sta]["components"].get(comp)
            if rec is None:
                continue
            y = top + ista * row_h + row_h / 2
            trace = raw_trace_on_grid(rec, start, end, dt)
            path = trace_path(normalize(trace), start, end, x0, y, panel_w, amp)
            parts.append(f'<path d="{path}" fill="none" stroke="#20242a" stroke-width="0.75"/>')
            if icomp == 0:
                parts.append(svg_text(x0 - 12, y + 4, station_axis_label(stations, sta, sort_by), 12, "end"))
        parts.append(f'<line x1="{x0:.1f}" y1="{height - bottom + 8:.1f}" x2="{x0 + panel_w:.1f}" y2="{height - bottom + 8:.1f}" stroke="#606975" stroke-width="1.2"/>')

    parts += [
        svg_text(width / 2, height - 12, "Time (s)", 16, color="#20242a"),
        "</svg>",
    ]
    return "\n".join(parts)


def make_grid(stations: dict, start: float, end: float, dt: float, sort_by: str) -> str:
    ordered = ordered_stations(stations, sort_by)
    width, height = 1500, 1850
    left, right, top, bottom = 130, 40, 105, 55
    col_gap = 26
    row_gap = 8
    col_w = (width - left - right - 2 * col_gap) / 3
    row_h = (height - top - bottom - (len(ordered) - 1) * row_gap) / len(ordered)
    amp = row_h * 0.34
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 35, "Raw Observed Displacement Waveforms", 24),
        svg_text(width / 2, 66, f"Unfiltered traces from disp/, sorted by {sort_by}", 13, color="#5f6872"),
    ]
    for icomp, comp in enumerate(COMPONENTS):
        x0 = left + icomp * (col_w + col_gap)
        parts.append(svg_text(x0 + col_w / 2, top - 15, COMP_LABELS.get(comp, comp), 15))
        for ista, sta in enumerate(ordered):
            y0 = top + ista * (row_h + row_gap)
            y = y0 + row_h / 2
            parts.append(f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{col_w:.1f}" height="{row_h:.1f}" fill="#f7f8fa" stroke="#e5e7eb"/>')
            rec = stations[sta]["components"].get(comp)
            if rec is not None:
                trace = raw_trace_on_grid(rec, start, end, dt)
                path = trace_path(normalize(trace), start, end, x0 + 4, y, col_w - 8, amp)
                parts.append(f'<path d="{path}" fill="none" stroke="#20242a" stroke-width="0.65"/>')
            if icomp == 0:
                parts.append(svg_text(x0 - 10, y + 4, sta, 10, "end"))
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--start", type=float, default=None, help="start time in seconds; default is the first SAC sample time")
    parser.add_argument("--end", type=float, default=None, help="end time in seconds; default is the last SAC sample time")
    parser.add_argument("--sort-by", choices=("dist", "az", "name"), default="dist")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    stations = station_records()
    default_start, default_end = default_time_range(stations)
    is_default_window = args.start is None and args.end is None
    start = default_start if args.start is None else args.start
    end = default_end if args.end is None else args.end
    suffix = args.sort_by
    window_suffix = time_suffix(start, end, is_default_window)
    section_path = RESULTS / f"raw_waveforms_record_section_{suffix}{window_suffix}.svg"
    grid_path = RESULTS / f"raw_waveforms_grid_{suffix}{window_suffix}.svg"
    section_path.write_text(make_record_section(stations, start, end, args.dt, args.sort_by))
    grid_path.write_text(make_grid(stations, start, end, args.dt, args.sort_by))
    print(f"Wrote {section_path}")
    print(f"Wrote {grid_path}")


if __name__ == "__main__":
    main()
