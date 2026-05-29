#!/usr/bin/env python3
"""Refine time shifts for the top grid-search models at sub-second spacing."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from array import array
from collections import defaultdict
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


def read_top_models(limit: int) -> list[dict]:
    with (RESULTS / "grid_search_results.csv").open() as f:
        rows = list(csv.DictReader(f))
    return [
        {
            "depth": int(row["depth_km"]),
            "strike": int(row["strike"]),
            "dip": int(row["dip"]),
            "rake": int(row["rake"]),
            "coarse_cc": float(row["mean_windowed_cc"]),
            "coarse_median_shift": float(row["median_shift_s"]),
        }
        for row in rows[:limit]
    ]


def coarse_shift_lookup() -> dict[tuple[int, int, int, int, str, str, str], float]:
    lookup = {}
    path = RESULTS / "trace_scores.csv"
    if not path.exists():
        return lookup
    with path.open() as f:
        for row in csv.DictReader(f):
            key = (
                int(row["depth_km"]),
                int(row["strike"]),
                int(row["dip"]),
                int(row["rake"]),
                row["station"],
                row["window"],
                row["component"],
            )
            lookup[key] = float(row["shift_s"])
    return lookup


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


def refine_corr(obs: array, syn: array, center_shift_s: float, half_width_s: float, shift_dt: float, sample_dt: float) -> tuple[float, float]:
    n_steps = int(round((2.0 * half_width_s) / shift_dt))
    start = center_shift_s - half_width_s
    best_cc = -2.0
    best_shift = center_shift_s
    seen = set()
    for i in range(n_steps + 1):
        shift_s = start + i * shift_dt
        shift_samples = int(round(shift_s / sample_dt))
        if shift_samples in seen:
            continue
        seen.add(shift_samples)
        cc = corr_at_shift(obs, syn, shift_samples)
        if cc is not None and cc > best_cc:
            best_cc = cc
            best_shift = shift_samples * sample_dt
    return best_cc, best_shift


def synthetic_components(sta: str, info: dict, model: dict, window: dict, sample_dt: float, band: tuple[float, float], mw: float, duration: float) -> dict[str, array]:
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
    sz = process_trace(syn["z"], window["start"], window["end"], sample_dt, band)
    sr = process_trace(syn["r"], window["start"], window["end"], sample_dt, band)
    st = process_trace(syn["t"], window["start"], window["end"], sample_dt, band)
    sn, se = rotate_rt_to_ne(sr, st, info["az"])
    return {"BHZ": sz, "BH1": sn, "BH2": se}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--sample-dt", type=float, default=0.1)
    parser.add_argument("--shift-dt", type=float, default=0.1)
    parser.add_argument("--half-width", type=float, default=1.0, help="fine search half-width around coarse shift in seconds")
    parser.add_argument("--bandpass", default="0.02/0.1")
    parser.add_argument("--mw", type=float, default=4.79)
    parser.add_argument("--duration", type=float, default=0.5)
    args = parser.parse_args()

    band = parse_band(args.bandpass)
    if band is None:
        raise SystemExit("bandpass must not be off")
    windows = [
        {"name": "pnl", "start": 0.0, "end": 90.0, "weight": 1.0},
        {"name": "surface", "start": 70.0, "end": 200.0, "weight": 1.0},
    ]
    stations = station_records()
    distances = [round(info["dist"], 1) for info in stations.values()]
    write_fk_model(MODEL_FILE, WORK / "Idaho_model")

    coarse = coarse_shift_lookup()
    details = []
    summary = []
    for rank, model in enumerate(read_top_models(args.top), 1):
        ensure_fk(model["depth"], distances, 2048, 0.1)
        ccs = []
        shifts = []
        for sta, info in stations.items():
            for window in windows:
                syn_by_comp = synthetic_components(sta, info, model, window, args.sample_dt, band, args.mw, args.duration)
                for comp in COMPONENTS:
                    obs = process_trace(info["components"][comp], window["start"], window["end"], args.sample_dt, band)
                    key = (model["depth"], model["strike"], model["dip"], model["rake"], sta, window["name"], comp)
                    center = coarse.get(key, model["coarse_median_shift"])
                    cc, shift = refine_corr(obs, syn_by_comp[comp], center, args.half_width, args.shift_dt, args.sample_dt)
                    ccs.append(cc)
                    shifts.append(shift)
                    details.append(
                        {
                            "rank": rank,
                            "depth_km": model["depth"],
                            "strike": model["strike"],
                            "dip": model["dip"],
                            "rake": model["rake"],
                            "station": sta,
                            "window": window["name"],
                            "component": comp,
                            "fine_cc": cc,
                            "fine_shift_s": shift,
                            "coarse_shift_s": center,
                        }
                    )
        score = sum(ccs) / len(ccs)
        med_shift = sorted(shifts)[len(shifts) // 2]
        summary.append(
            {
                "rank": rank,
                "depth_km": model["depth"],
                "strike": model["strike"],
                "dip": model["dip"],
                "rake": model["rake"],
                "coarse_mean_windowed_cc": model["coarse_cc"],
                "fine_mean_windowed_cc": score,
                "fine_median_shift_s": med_shift,
                "n_window_traces": len(ccs),
            }
        )
        print(
            f"rank {rank}: fine cc={score:.4f} depth={model['depth']} "
            f"strike={model['strike']} dip={model['dip']} rake={model['rake']}"
        )

    summary_path = RESULTS / "refined_top_model_scores.csv"
    details_path = RESULTS / "refined_trace_scores_top5.csv"
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    with details_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(details[0].keys()))
        writer.writeheader()
        writer.writerows(details)
    print(f"Wrote {summary_path}")
    print(f"Wrote {details_path}")


if __name__ == "__main__":
    main()
