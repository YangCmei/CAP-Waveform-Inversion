#!/usr/bin/env python3
"""Plot filtered observed displacement waveforms for presentation figures."""

from __future__ import annotations

import argparse
import math
import sys
from array import array
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from waveform_inversion import COMPONENTS, RESULTS, parse_band, process_trace, station_records


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


def normalize(values: array) -> list[float]:
    if not values:
        return []
    scale = max(abs(v) for v in values)
    if scale <= 0.0 or not math.isfinite(scale):
        return [0.0 for _ in values]
    return [v / scale for v in values]


def trace_path(values: list[float], start: float, end: float, x0: float, y0: float, width: float, amp: float, max_points: int = 900) -> str:
    if not values:
        return ""
    n = len(values)
    step = max(1, n // max_points)
    denom = max(1, n - 1)
    pts = []
    for i in range(0, n, step):
        t = start + (end - start) * i / denom
        x = x0 + (t - start) / (end - start) * width
        y = y0 - values[i] * amp
        pts.append((x, y))
    if pts[-1][0] < x0 + width:
        pts.append((x0 + width, y0 - values[-1] * amp))
    return "M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts)


def station_order(stations: dict) -> list[str]:
    return sorted(stations, key=lambda name: stations[name]["dist"])


def default_end_time(stations: dict) -> float:
    return max(
        rec["e"]
        for info in stations.values()
        for rec in info["components"].values()
    )


def time_suffix(start: float, end: float, is_default_end: bool) -> str:
    if is_default_end and start == 0.0:
        return ""
    return f"_{start:g}_{end:g}s".replace(".", "p")


def make_record_section(stations: dict, start: float, end: float, dt: float, band: tuple[float, float]) -> str:
    ordered = station_order(stations)
    width, height = 1500, 930
    margin_left, margin_right = 180, 40
    top, bottom = 135, 65
    panel_gap = 36
    panel_w = (width - margin_left - margin_right - 2 * panel_gap) / 3
    row_h = (height - top - bottom) / len(ordered)
    amp = row_h * 0.36

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 35, "Filtered Observed Waveforms Used in Scoring", 24),
        svg_text(width / 2, 66, f"Band-pass {band[0]:g}-{band[1]:g} Hz; resampled to {dt:g} s; each trace normalized independently", 13, color="#5f6872"),
    ]

    for icomp, comp in enumerate(COMPONENTS):
        x0 = margin_left + icomp * (panel_w + panel_gap)
        parts += [
            f'<rect x="{x0:.1f}" y="{top - 24:.1f}" width="{panel_w:.1f}" height="{height - top - bottom + 28:.1f}" fill="#f7f8fa" stroke="#d8dde3"/>',
            svg_text(x0 + panel_w / 2, top - 38, COMP_LABELS.get(comp, comp), 15),
            f'<line x1="{x0:.1f}" y1="{height - bottom + 8:.1f}" x2="{x0 + panel_w:.1f}" y2="{height - bottom + 8:.1f}" stroke="#9aa3ad"/>',
            svg_text(x0, height - bottom + 24, f"{start:g}s", 11, "start", "#5f6872"),
            svg_text(x0 + panel_w, height - bottom + 24, f"{end:g}s", 11, "end", "#5f6872"),
        ]
        for tick in range(int(start), int(end) + 1, 50):
            x = x0 + (tick - start) / (end - start) * panel_w
            parts.append(f'<line x1="{x:.1f}" y1="{top - 12:.1f}" x2="{x:.1f}" y2="{height - bottom + 8:.1f}" stroke="#e5e7eb" stroke-width="1"/>')

        for ista, sta in enumerate(ordered):
            rec = stations[sta]["components"].get(comp)
            if rec is None:
                continue
            y = top + ista * row_h + row_h / 2
            trace = process_trace(rec, start, end, dt, band)
            path = trace_path(normalize(trace), start, end, x0, y, panel_w, amp)
            parts.append(f'<path d="{path}" fill="none" stroke="#20242a" stroke-width="0.9"/>')
            if icomp == 0:
                parts.append(svg_text(x0 - 12, y + 4, f'{sta} ({stations[sta]["dist"]:.0f} km)', 10, "end", "#20242a"))

    parts += [
        svg_text(width / 2, height - 16, "Time (s)", 12, color="#5f6872"),
        "</svg>",
    ]
    return "\n".join(parts)


def make_raw_vs_filtered_examples(stations: dict, start: float, end: float, dt: float, band: tuple[float, float]) -> str:
    ordered = station_order(stations)
    picks = [ordered[0], ordered[len(ordered) // 2], ordered[-1]]
    width, height = 1500, 720
    left, right, top = 170, 40, 95
    panel_gap = 36
    panel_w = (width - left - right - 2 * panel_gap) / 3
    row_h = 170
    amp = 48
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 35, "Low-frequency filtered waveforms", 24),
        svg_text(width / 2, 62, f"Gray: raw displacement, blue: {band[0]:g}-{band[1]:g} Hz filtered; each trace normalized independently", 13, color="#5f6872"),
    ]
    for icomp, comp in enumerate(COMPONENTS):
        x0 = left + icomp * (panel_w + panel_gap)
        parts.append(svg_text(x0 + panel_w / 2, 88, COMP_LABELS.get(comp, comp), 15))
        for ipick, sta in enumerate(picks):
            rec = stations[sta]["components"].get(comp)
            if rec is None:
                continue
            y = top + ipick * row_h + 70
            raw = process_trace(rec, start, end, dt, None)
            filt = process_trace(rec, start, end, dt, band)
            raw_path = trace_path(normalize(raw), start, end, x0, y, panel_w, amp)
            filt_path = trace_path(normalize(filt), start, end, x0, y, panel_w, amp)
            parts += [
                f'<rect x="{x0:.1f}" y="{y - 62:.1f}" width="{panel_w:.1f}" height="124" fill="#f7f8fa" stroke="#d8dde3"/>',
                f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x0 + panel_w:.1f}" y2="{y:.1f}" stroke="#d8dde3"/>',
                f'<path d="{raw_path}" fill="none" stroke="#a4acb5" stroke-width="0.8"/>',
                f'<path d="{filt_path}" fill="none" stroke="#1f6f8b" stroke-width="1.3"/>',
            ]
            if icomp == 0:
                parts.append(svg_text(x0 - 12, y + 4, f'{sta} ({stations[sta]["dist"]:.0f} km)', 11, "end"))
        time_y = top + 2 * row_h + 70 + 84
        parts.append(svg_text(x0, time_y, f"{start:g}s", 11, "start", "#5f6872"))
        parts.append(svg_text(x0 + panel_w, time_y, f"{end:g}s", 11, "end", "#5f6872"))
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bandpass", default="0.02/0.1", help="band-pass in Hz as fmin/fmax")
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--start", type=float, default=0.0)
    parser.add_argument("--end", type=float, default=None, help="end time in seconds; default is the full record length")
    args = parser.parse_args()

    band = parse_band(args.bandpass)
    if band is None:
        raise SystemExit("--bandpass must not be off for this plotting script")
    RESULTS.mkdir(exist_ok=True)
    stations = station_records()
    is_default_end = args.end is None
    end = default_end_time(stations) if is_default_end else args.end
    record = make_record_section(stations, args.start, end, args.dt, band)
    examples = make_raw_vs_filtered_examples(stations, args.start, end, args.dt, band)
    suffix = time_suffix(args.start, end, is_default_end)
    record_path = RESULTS / f"filtered_waveforms_record_section{suffix}.svg"
    examples_path = RESULTS / f"low_frequency_waveform_examples{suffix}.svg"
    record_path.write_text(record)
    examples_path.write_text(examples)
    print(f"Wrote {record_path}")
    print(f"Wrote {examples_path}")


if __name__ == "__main__":
    main()
