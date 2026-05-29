# Waveform Inversion Result

## Data and Method

- Observed data: 12 stations x 3 components, instrument-response-removed displacement SAC files in `data/disp/`.
- Synthetics: FK Green functions plus `syn`, using the supplied Idaho 1-D crustal model.
- Moment magnitude: fixed at `Mw = 4.79`.
- Source time function: fixed as `syn -D0.5`, a simple 0.5 s source-time function used by the FK `syn` program.
- Band-pass: both observed and synthetic waveforms are filtered to `0.02-0.1 Hz` before scoring.
- Windows: simplified CAP-like `Pnl: 0-90 s` and `surface: 70-200 s`.
- Time shift is not a source-model parameter. It is only an auxiliary alignment value.
- In the corrected re-score, each station-window tests shifts from `2` to `6 s` at `0.1 s` spacing. The same shift is used for all 3 components of that station-window. Different stations and different windows may have different shifts.
- Objective function: for each station-window, choose the shift that maximizes the mean normalized correlation of the 3 components, then average the `12 stations x 2 windows = 24` station-window scores. This gives `Mean_windowed_CC`.

If written as an error-like quantity:

```text
Misfit = 1 - Mean_windowed_CC
```

## Search Workflow

The initial reference solution in `models/evt_info.txt` is:

```text
Depth  = 12 km
Strike = 330 deg
Dip    = 80 deg
Rake   = -15 deg
Mw     = 4.79
```

The first coarse grid was:

```text
Depth  = 6, 10, 14, 18 km
Strike = 0, 60, 120, 180, 240, 300 deg
Dip    = 30, 60, 90 deg
Rake   = -180, -120, -60, 0, 60, 120, 180 deg
```

The focused refinement around the coarse high-score region was:

```text
Depth  = 4, 6, 8, 10 km
Strike = 40, 50, 60, 70, 80, 90 deg
Dip    = 40, 50, 60, 70, 80, 90 deg
Rake   = -180, -170, -160, -150, -140, -130, -120, -110, -100, -90 deg
```

The original focused grid has `1440` models and is stored in `results/grid_search_results.csv`. After noticing that the `1 s` shift spacing was too coarse, the original top 20 candidates were re-scored with `2-6 s` shifts at `0.1 s` spacing using the corrected station-window shared-shift rule.

Strictly speaking, the current corrected result is the best model within that re-scored top-20 candidate set, not a full `0.1 s` re-score of all 1440 models.

## Best Model

Current best corrected candidate:

```text
Depth  = 6 km
Strike = 70 deg
Dip    = 60 deg
Rake   = -160 deg
Mean_windowed_CC = 0.6844
```

The score is averaged over 24 station-window scores. Within each station-window, the 3 components share one auxiliary shift and their correlations are averaged.

## Top 10 Corrected Candidates

| Rank | Depth (km) | Strike | Dip | Rake | Mean_windowed_CC | n station-window scores |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 6 | 70 | 60 | -160 | 0.6844 | 24 |
| 2 | 6 | 70 | 50 | -160 | 0.6727 | 24 |
| 3 | 6 | 70 | 60 | -170 | 0.6577 | 24 |
| 4 | 6 | 70 | 50 | -170 | 0.6449 | 24 |
| 5 | 6 | 60 | 70 | -160 | 0.6435 | 24 |
| 6 | 6 | 60 | 60 | -160 | 0.6423 | 24 |
| 7 | 6 | 60 | 70 | -150 | 0.6359 | 24 |
| 8 | 6 | 60 | 70 | -170 | 0.6358 | 24 |
| 9 | 6 | 60 | 60 | -170 | 0.6329 | 24 |
| 10 | 6 | 60 | 80 | -160 | 0.6264 | 24 |

The more stable interpretation is a high-score region:

```text
Depth  around 6 km
Strike around 60-70 deg
Dip    around 50-70 deg
Rake   around -170 to -150 deg
```

## Interpretation

The corrected result is closer to the USGS first nodal plane (`66/78/-148`) than the previous `1 s`-shift result. The strike is very similar, and the rake now falls in a nearby mostly strike-slip range. The remaining depth difference is larger: this simplified workflow prefers a shallow `6 km` candidate, while USGS reports about `12.2 km`.

That depth difference should be interpreted cautiously. The workflow uses a 1-D velocity model, fixed time windows, no station/path corrections, no station or azimuthal weighting, and only the original top 20 candidates were re-scored at `0.1 s`. The current result is therefore best stated as a preferred high-score region under this simplified scoring setup, not as a unique global solution.

## Output Files

- `results/grid_search_results.csv`: original 1440-model focused grid with the older `1 s` shift scoring.
- `results/trace_scores.csv`: trace-level details for the older scoring.
- `results/grid_search_results_shift2_6_dt0p1.csv`: corrected `2-6 s`, `0.1 s` re-score of the original top 20 candidates.
- `results/trace_scores_shift2_6_dt0p1.csv`: station/window/component CC details for the corrected re-score. Components in the same station-window share the same auxiliary shift.
- `results/waveform_fit_top5_index.md`: index of the top-5 observed/synthetic waveform comparison figures.
- `results/shift2_6_best_model_station_scores.svg`: station-score map for the corrected best candidate.

## Reproducibility

Coarse grid:

```bash
python3 scripts/waveform_inversion.py --depth-min 6 --depth-max 18 --depth-step 4 --strike-min 0 --strike-max 300 --strike-step 60 --dip-min 30 --dip-max 90 --dip-step 30 --rake-min -180 --rake-max 180 --rake-step 60 --score-dt 1.0 --max-shift 40 --windows pnl:0:90:1,surface:70:200:1 --bandpass 0.02/0.1
```

Focused refinement:

```bash
python3 scripts/waveform_inversion.py --depth-min 4 --depth-max 10 --depth-step 2 --strike-min 40 --strike-max 90 --strike-step 10 --dip-min 40 --dip-max 90 --dip-step 10 --rake-min -180 --rake-max -90 --rake-step 10 --score-dt 1.0 --max-shift 40 --windows pnl:0:90:1,surface:70:200:1 --bandpass 0.02/0.1 --keep-details
```

Corrected candidate re-score:

```bash
python3 scripts/rescore_grid_fine_shift.py --sample-dt 0.1 --shift-min 2 --shift-max 6 --shift-dt 0.1 --bandpass 0.02/0.1 --limit 20
```
