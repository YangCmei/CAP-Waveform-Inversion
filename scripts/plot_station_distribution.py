#!/usr/bin/env python3
"""Plot station distribution and azimuthal coverage without extra packages."""

from __future__ import annotations

import csv
import math
import struct
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


def station_metadata() -> list[dict]:
    stations = {}
    for path in sorted(DISP_DIR.glob("*.sac.dis")):
        header = read_sac_header(path)
        if header["component"] != "BHZ":
            continue
        network, station = path.name.split(".")[:2]
        key = f"{network}.{station}"
        stations[key] = header | {"name": key}
    return list(stations.values())


def polar_to_xy(distance_km: float, az_deg: float) -> tuple[float, float]:
    az = math.radians(az_deg)
    return distance_km * math.sin(az), distance_km * math.cos(az)


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 13,
    anchor: str = "middle",
    color: str = "#20242a",
    weight: str = "400",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Microsoft YaHei, sans-serif" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{color}">{text}</text>'
    )


def make_svg(stations: list[dict]) -> str:
    width, height = 1400, 980
    left_cx, left_cy = 360, 390
    right_cx, right_cy = 1040, 390
    max_dist = max(s["dist"] for s in stations)
    ring_step = 100
    ring_max = int(math.ceil(max_dist / ring_step) * ring_step)
    radius = 250
    scale = radius / ring_max

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f7f8fa"/>',
        svg_text(width / 2, 42, "Station Distribution and Azimuthal Coverage", 26, weight="700"),
        svg_text(left_cx, 88, "Map view relative to epicenter", 17, weight="700"),
        svg_text(right_cx, 88, "Polar view: azimuth and epicentral distance", 17, weight="700"),
    ]

    for cx, cy in ((left_cx, left_cy), (right_cx, right_cy)):
        for dist in range(ring_step, ring_max + 1, ring_step):
            r = dist * scale
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r:.1f}" fill="none" stroke="#d8dde3" stroke-width="1"/>')
            parts.append(svg_text(cx + 5, cy - r + 14, f"{dist} km", 10, "start"))
        parts.append(f'<line x1="{cx-radius}" y1="{cy}" x2="{cx+radius}" y2="{cy}" stroke="#c2c8cf" stroke-width="1"/>')
        parts.append(f'<line x1="{cx}" y1="{cy-radius}" x2="{cx}" y2="{cy+radius}" stroke="#c2c8cf" stroke-width="1"/>')
        parts.append(svg_text(cx, cy - radius - 18, "N", 15))
        parts.append(svg_text(cx + radius + 18, cy + 5, "E", 15, "start"))
        parts.append(svg_text(cx, cy + radius + 28, "S", 15))
        parts.append(svg_text(cx - radius - 18, cy + 5, "W", 15, "end"))
        parts.append(f'<circle cx="{cx}" cy="{cy}" r="6" fill="#a64b2a"/>')
        parts.append(svg_text(cx, cy - 12, "Epicenter", 11))

    for az in range(0, 360, 45):
        east, north = polar_to_xy(radius, az)
        x2 = right_cx + east
        y2 = right_cy - north
        parts.append(
            f'<line x1="{right_cx:.1f}" y1="{right_cy:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            'stroke="#b8c0c8" stroke-width="0.8" stroke-dasharray="4 4"/>'
        )
        lx = right_cx + (radius + 42) * math.sin(math.radians(az))
        ly = right_cy - (radius + 42) * math.cos(math.radians(az)) + 5
        anchor = "middle"
        if 0 < az < 180:
            anchor = "start"
        elif 180 < az < 360:
            anchor = "end"
        parts.append(svg_text(lx, ly, f"{az}°", 15, anchor, "#20242a", "700"))

    colors = ["#1f6f8b", "#5d7a3b", "#8a6f2a", "#595e73", "#a64b2a"]
    for i, sta in enumerate(stations, 1):
        east, north = polar_to_xy(sta["dist"], sta["az"])
        x = left_cx + east * scale
        y = left_cy - north * scale
        color = colors[(i - 1) % len(colors)]
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" fill="{color}" stroke="#ffffff" stroke-width="1.5"/>')
        parts.append(svg_text(x, y + 3.5, str(i), 9, "middle", "#ffffff", "700"))

        # Polar panel uses the same projection from azimuth/distance; draw a ray
        px = right_cx + east * scale
        py = right_cy - north * scale
        parts.append(f'<line x1="{right_cx}" y1="{right_cy}" x2="{px:.1f}" y2="{py:.1f}" stroke="{color}" stroke-width="1" opacity="0.35"/>')
        parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="8" fill="{color}" stroke="#ffffff" stroke-width="1.5"/>')
        parts.append(svg_text(px, py + 3.5, str(i), 9, "middle", "#ffffff", "700"))

    az_sorted = sorted(s["az"] for s in stations)
    gaps = []
    for i, az in enumerate(az_sorted):
        nxt = az_sorted[(i + 1) % len(az_sorted)]
        if i == len(az_sorted) - 1:
            nxt += 360.0
        gaps.append(nxt - az)
    max_gap = max(gaps)
    min_dist = min(s["dist"] for s in stations)
    parts += [
        '<rect x="80" y="730" width="1240" height="58" rx="4" fill="#ffffff" stroke="#d8dde3"/>',
        svg_text(110, 766, f"Stations: {len(stations)}", 16, "start", weight="700"),
        svg_text(300, 766, f"Distance range: {min_dist:.1f}-{max_dist:.1f} km", 16, "start", weight="700"),
        svg_text(650, 766, f"Largest azimuthal gap: {max_gap:.1f}°", 16, "start", weight="700"),
        svg_text(995, 766, "Smaller gap means better coverage", 16, "start", weight="700"),
        svg_text(110, 835, "Station legend", 17, "start", weight="700"),
    ]
    col_x = [110, 415, 720, 1025]
    row_y = [865, 896, 927]
    for i, sta in enumerate(stations, 1):
        col = (i - 1) // 3
        row = (i - 1) % 3
        text = f'{i}. {sta["name"]}  d={sta["dist"]:.0f} km  az={sta["az"]:.0f}°'
        parts.append(svg_text(col_x[col], row_y[row], text, 15, "start", "#111827", "700"))
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    stations = station_metadata()
    with (RESULTS / "station_distribution.csv").open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["name", "station", "stla", "stlo", "evla", "evlo", "dist", "az", "baz"],
        )
        writer.writeheader()
        for sta in stations:
            writer.writerow({key: sta[key] for key in writer.fieldnames})
    (RESULTS / "station_distribution.svg").write_text(make_svg(stations))
    print(f"Wrote {RESULTS / 'station_distribution.svg'}")
    print(f"Wrote {RESULTS / 'station_distribution.csv'}")


if __name__ == "__main__":
    main()
