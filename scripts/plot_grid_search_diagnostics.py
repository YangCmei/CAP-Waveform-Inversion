#!/usr/bin/env python3
"""Plot diagnostic figures from grid_search_results.csv."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def svg_text(x, y, text, size=12, anchor="middle", color="#20242a", weight="400"):
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Microsoft YaHei, sans-serif" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{color}">{text}</text>'
    )


def read_rows(path: Path):
    with path.open() as f:
        rows = []
        for row in csv.DictReader(f):
            rows.append(
                {
                    "depth": int(row["depth_km"]),
                    "strike": int(row["strike"]),
                    "dip": int(row["dip"]),
                    "rake": int(row["rake"]),
                    "cc": float(row["mean_windowed_cc"]),
                    "shift": float(row["median_shift_s"]),
                }
            )
    return rows


def color_scale(value, vmin, vmax):
    if vmax <= vmin:
        t = 1.0
    else:
        t = (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    # light gray-blue to orange-red
    c0 = (232, 243, 248)
    c1 = (166, 75, 42)
    rgb = tuple(round(c0[i] * (1 - t) + c1[i] * t) for i in range(3))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def aggregate_best(rows, x_key, y_key):
    best = {}
    for row in rows:
        key = (row[x_key], row[y_key])
        if key not in best or row["cc"] > best[key]["cc"]:
            best[key] = row
    xs = sorted({row[x_key] for row in rows})
    ys = sorted({row[y_key] for row in rows})
    return xs, ys, best


def heatmap(rows, x_key, y_key, title, x_label, y_label, out_path):
    xs, ys, best = aggregate_best(rows, x_key, y_key)
    values = [row["cc"] for row in best.values()]
    vmin, vmax = min(values), max(values)
    width, height = 1050, 760
    left, right, top, bottom = 165, 145, 85, 105
    grid_w = width - left - right
    grid_h = height - top - bottom
    cell_w = grid_w / len(xs)
    cell_h = grid_h / len(ys)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 38, title, 25, weight="700"),
    ]
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            row = best.get((x, y))
            xx = left + ix * cell_w
            yy = top + iy * cell_h
            fill = color_scale(row["cc"], vmin, vmax) if row else "#f1f3f5"
            parts.append(f'<rect x="{xx:.1f}" y="{yy:.1f}" width="{cell_w:.1f}" height="{cell_h:.1f}" fill="{fill}" stroke="#ffffff" stroke-width="1"/>')
            if row:
                parts.append(svg_text(xx + cell_w / 2, yy + cell_h / 2 + 5, f'{row["cc"]:.3f}', 12, color="#111827", weight="700"))
    for ix, x in enumerate(xs):
        parts.append(svg_text(left + ix * cell_w + cell_w / 2, top + grid_h + 28, str(x), 13, weight="700"))
    for iy, y in enumerate(ys):
        parts.append(svg_text(left - 15, top + iy * cell_h + cell_h / 2 + 5, str(y), 13, "end", weight="700"))
    parts += [
        f'<rect x="{left:.1f}" y="{top:.1f}" width="{grid_w:.1f}" height="{grid_h:.1f}" fill="none" stroke="#20242a" stroke-width="1.2"/>',
        svg_text(left + grid_w / 2, height - 28, x_label, 17, weight="700"),
        svg_text(65, top + grid_h / 2, y_label, 17, weight="700"),
    ]

    # color legend
    lx, ly, lw, lh = width - 100, top, 24, grid_h
    for i in range(80):
        t0 = i / 80
        val = vmax - t0 * (vmax - vmin)
        yy = ly + i * lh / 80
        parts.append(f'<rect x="{lx}" y="{yy:.1f}" width="{lw}" height="{lh/80 + 1:.1f}" fill="{color_scale(val, vmin, vmax)}" stroke="none"/>')
    parts.append(f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" fill="none" stroke="#6b7280"/>')
    parts.append(svg_text(lx + lw + 10, ly + 5, f"{vmax:.3f}", 12, "start", weight="700"))
    parts.append(svg_text(lx + lw + 10, ly + lh, f"{vmin:.3f}", 12, "start", weight="700"))
    parts.append("</svg>")
    out_path.write_text("\n".join(parts))


def top_bar(rows, out_path, n=20):
    rows = sorted(rows, key=lambda r: r["cc"], reverse=True)[:n]
    width, height = 1300, 760
    left, right, top, bottom = 125, 40, 105, 185
    plot_w = width - left - right
    plot_h = height - top - bottom
    vmin = min(r["cc"] for r in rows) - 0.002
    vmax = max(r["cc"] for r in rows)
    bar_gap = 6
    bar_w = (plot_w - bar_gap * (n - 1)) / n
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 38, f"Top {n} Grid-search Models", 25, weight="700"),
        svg_text(left, 78, "Mean_windowed_CC", 17, "start", weight="700"),
    ]
    for i, row in enumerate(rows):
        x = left + i * (bar_w + bar_gap)
        h = (row["cc"] - vmin) / (vmax - vmin) * plot_h
        y = top + plot_h - h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="#1f6f8b"/>')
        parts.append(svg_text(x + bar_w / 2, y - 6, f'{row["cc"]:.3f}', 10, weight="700"))
        label = f'{row["depth"]}/{row["strike"]}/{row["dip"]}/{row["rake"]}'
        # rotate label
        tx = x + bar_w / 2
        ty = top + plot_h + 18
        parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="10" font-family="Arial, Microsoft YaHei, sans-serif" '
            f'text-anchor="end" font-weight="700" fill="#20242a" transform="rotate(-55 {tx:.1f} {ty:.1f})">{label}</text>'
        )
    parts += [
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#20242a" stroke-width="1.2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#20242a" stroke-width="1.2"/>',
        svg_text(width / 2, height - 24, "Depth / Strike / Dip / Rake", 15, weight="700"),
        "</svg>",
    ]
    out_path.write_text("\n".join(parts))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=str(RESULTS / "grid_search_results.csv"))
    args = parser.parse_args()
    rows = read_rows(Path(args.results))
    heatmap(rows, "strike", "depth", "Best Mean_windowed_CC by Depth and Strike", "Strike (deg)", "Depth (km)", RESULTS / "heatmap_depth_strike.svg")
    heatmap(rows, "rake", "dip", "Best Mean_windowed_CC by Dip and Rake", "Rake (deg)", "Dip (deg)", RESULTS / "heatmap_dip_rake.svg")
    heatmap(rows, "rake", "strike", "Best Mean_windowed_CC by Strike and Rake", "Rake (deg)", "Strike (deg)", RESULTS / "heatmap_strike_rake.svg")
    top_bar(rows, RESULTS / "top20_grid_models.svg")
    print(f"Wrote {RESULTS / 'heatmap_depth_strike.svg'}")
    print(f"Wrote {RESULTS / 'heatmap_dip_rake.svg'}")
    print(f"Wrote {RESULTS / 'heatmap_strike_rake.svg'}")
    print(f"Wrote {RESULTS / 'top20_grid_models.svg'}")


if __name__ == "__main__":
    main()
