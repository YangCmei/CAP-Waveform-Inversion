#!/usr/bin/env python3
"""Plot station-wise scores for a selected focal-mechanism model."""

from __future__ import annotations

import argparse
import csv
import math
import struct
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DISP_DIR = next(path for path in (ROOT / "data" / "disp", ROOT / "disp") if path.exists())
RESULTS = ROOT / "results"


def read_sac_header(path: Path) -> dict:
    raw = path.read_bytes()[:632]
    floats = struct.unpack("<70f", raw[:280])
    chars = raw[440:632]
    keys = [
        chars[i * 8 : (i + 1) * 8].decode("ascii", "ignore").strip("\x00 ")
        for i in range(24)
    ]
    return {
        "station": keys[0],
        "component": keys[20],
        "evla": floats[35],
        "evlo": floats[36],
        "stla": floats[31],
        "stlo": floats[32],
        "dist": floats[50],
        "az": floats[51],
        "baz": floats[52],
    }


def station_metadata() -> dict[str, dict]:
    stations = {}
    for path in sorted(DISP_DIR.glob("*.sac.dis")):
        header = read_sac_header(path)
        if header["component"] != "BHZ":
            continue
        network, station = path.name.split(".")[:2]
        key = f"{network}.{station}"
        stations[key] = header | {"name": key}
    return stations


def station_scores(depth: int, strike: int, dip: int, rake: int) -> dict[str, dict]:
    fine = scores_from_file(RESULTS / "trace_scores_shift2_6_dt0p1.csv", depth, strike, dip, rake, "cc", "shift_s")
    if fine:
        return fine
    refined = scores_from_file(RESULTS / "refined_trace_scores_top5.csv", depth, strike, dip, rake, "fine_cc", "fine_shift_s")
    if refined:
        return refined
    return scores_from_file(RESULTS / "trace_scores.csv", depth, strike, dip, rake, "cc", "shift_s")


def scores_from_file(path: Path, depth: int, strike: int, dip: int, rake: int, cc_key: str, shift_key: str) -> dict[str, dict]:
    if not path.exists():
        return {}
    scores: dict[str, list[float]] = defaultdict(list)
    shifts: dict[str, list[float]] = defaultdict(list)
    with path.open() as f:
        for row in csv.DictReader(f):
            if (
                int(row["depth_km"]) == depth
                and int(row["strike"]) == strike
                and int(row["dip"]) == dip
                and int(row["rake"]) == rake
            ):
                scores[row["station"]].append(float(row[cc_key]))
                shifts[row["station"]].append(float(row[shift_key]))
    return {
        sta: {
            "mean_cc": sum(vals) / len(vals),
            "min_cc": min(vals),
            "max_cc": max(vals),
            "n": len(vals),
            "median_shift": sorted(shifts[sta])[len(shifts[sta]) // 2],
        }
        for sta, vals in scores.items()
    }


def polar_to_xy(distance_km: float, az_deg: float) -> tuple[float, float]:
    az = math.radians(az_deg)
    return distance_km * math.sin(az), distance_km * math.cos(az)


def svg_text(x, y, text, size=12, anchor="middle", color="#20242a", weight="400"):
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Microsoft YaHei, sans-serif" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{color}">{text}</text>'
    )


def color_scale(value: float, vmin: float, vmax: float) -> str:
    if vmax <= vmin:
        t = 1.0
    else:
        t = (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    # lower score: light yellow; higher score: blue-green
    c0 = (255, 247, 230)
    c1 = (31, 111, 139)
    rgb = tuple(round(c0[i] * (1 - t) + c1[i] * t) for i in range(3))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


LABEL_OFFSETS = {
    "IW.PLID": (10, -16, "start"),
    "US.BMO": (-12, -16, "end"),
    "IW.MFID": (-12, 18, "end"),
    "IW.DLMT": (12, -18, "start"),
    "US.WVOR": (-12, 22, "end"),
    "US.MSO": (12, -20, "start"),
    "US.HAWA": (-12, -18, "end"),
    "US.BOZ": (12, 18, "start"),
    "US.ELK": (-12, 20, "end"),
    "US.NEW": (-12, -18, "end"),
    "US.HWUT": (12, -18, "start"),
    "US.DUG": (12, 20, "start"),
}


def make_svg(stations: dict[str, dict], scores: dict[str, dict], model: tuple[int, int, int, int]) -> str:
    width, height = 1120, 920
    cx, cy = 500, 430
    radius = 315
    max_dist = max(sta["dist"] for sta in stations.values())
    ring_step = 100
    ring_max = int(math.ceil(max_dist / ring_step) * ring_step)
    scale = radius / ring_max
    vmin = min(s["mean_cc"] for s in scores.values())
    vmax = max(s["mean_cc"] for s in scores.values())
    depth, strike, dip, rake = model

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(width / 2, 38, "Station Mean_windowed_CC for Best Model", 25, weight="700"),
        svg_text(width / 2, 68, f"Depth={depth} km, Strike={strike}°, Dip={dip}°, Rake={rake}°", 15, color="#5f6872", weight="700"),
    ]

    for dist in range(ring_step, ring_max + 1, ring_step):
        r = dist * scale
        parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r:.1f}" fill="none" stroke="#d8dde3" stroke-width="1"/>')
        parts.append(svg_text(cx + 5, cy - r + 14, f"{dist} km", 11, "start", "#5f6872"))
    parts += [
        f'<line x1="{cx-radius}" y1="{cy}" x2="{cx+radius}" y2="{cy}" stroke="#c2c8cf" stroke-width="1"/>',
        f'<line x1="{cx}" y1="{cy-radius}" x2="{cx}" y2="{cy+radius}" stroke="#c2c8cf" stroke-width="1"/>',
        svg_text(cx, cy - radius - 18, "N", 16, weight="700"),
        svg_text(cx + radius + 20, cy + 5, "E", 16, "start", weight="700"),
        svg_text(cx, cy + radius + 30, "S", 16, weight="700"),
        svg_text(cx - radius - 20, cy + 5, "W", 16, "end", weight="700"),
        f'<circle cx="{cx}" cy="{cy}" r="7" fill="#a64b2a"/>',
        svg_text(cx, cy - 14, "Epicenter", 12, color="#a64b2a", weight="700"),
    ]

    for sta, meta in sorted(stations.items(), key=lambda item: item[1]["az"]):
        if sta not in scores:
            continue
        east, north = polar_to_xy(meta["dist"], meta["az"])
        x = cx + east * scale
        y = cy - north * scale
        mean_cc = scores[sta]["mean_cc"]
        fill = color_scale(mean_cc, vmin, vmax)
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="12" fill="{fill}" stroke="#20242a" stroke-width="1.2"/>')
        dx, dy, anchor = LABEL_OFFSETS.get(sta, (12, -14, "start"))
        label = f"{sta}  {mean_cc:.2f}"
        parts.append(svg_text(x + dx, y + dy, label, 14, anchor, "#111827", "700"))

    # colorbar
    lx, ly, lw, lh = 910, 205, 28, 390
    parts.append(svg_text(lx + 15, ly - 22, "Mean CC", 15, weight="700"))
    for i in range(100):
        val = vmax - (i / 100) * (vmax - vmin)
        yy = ly + i * lh / 100
        parts.append(f'<rect x="{lx}" y="{yy:.1f}" width="{lw}" height="{lh/100 + 1:.1f}" fill="{color_scale(val, vmin, vmax)}" stroke="none"/>')
    parts.append(f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" fill="none" stroke="#6b7280"/>')
    parts.append(svg_text(lx + lw + 12, ly + 5, f"{vmax:.2f}", 13, "start", weight="700"))
    parts.append(svg_text(lx + lw + 12, ly + lh, f"{vmin:.2f}", 13, "start", weight="700"))
    parts.append("</svg>")
    return "\n".join(parts)


def write_csv(stations: dict[str, dict], scores: dict[str, dict], model: tuple[int, int, int, int], prefix: str) -> Path:
    out = RESULTS / f"{prefix}.csv"
    with out.open("w", newline="") as f:
        fields = ["station", "dist_km", "az_deg", "mean_cc", "min_cc", "max_cc", "n_scores", "median_shift_s", "depth", "strike", "dip", "rake"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for sta, meta in sorted(stations.items(), key=lambda item: item[1]["dist"]):
            s = scores.get(sta)
            if not s:
                continue
            writer.writerow(
                {
                    "station": sta,
                    "dist_km": f'{meta["dist"]:.1f}',
                    "az_deg": f'{meta["az"]:.1f}',
                    "mean_cc": f'{s["mean_cc"]:.6f}',
                    "min_cc": f'{s["min_cc"]:.6f}',
                    "max_cc": f'{s["max_cc"]:.6f}',
                    "n_scores": s["n"],
                    "median_shift_s": f'{s["median_shift"]:.1f}',
                    "depth": model[0],
                    "strike": model[1],
                    "dip": model[2],
                    "rake": model[3],
                }
            )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", type=int, default=6)
    parser.add_argument("--strike", type=int, default=60)
    parser.add_argument("--dip", type=int, default=70)
    parser.add_argument("--rake", type=int, default=-180)
    parser.add_argument("--output-prefix", default="best_model_station_scores")
    args = parser.parse_args()
    model = (args.depth, args.strike, args.dip, args.rake)
    stations = station_metadata()
    scores = station_scores(*model)
    if not scores:
        raise SystemExit(f"no scores found for model {model}")
    svg = make_svg(stations, scores, model)
    out_svg = RESULTS / f"{args.output_prefix}.svg"
    out_svg.write_text(svg)
    out_csv = write_csv(stations, scores, model, args.output_prefix)
    print(f"Wrote {out_svg}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
