# 02_preprocess_eeglab — EDF → preprocessing → PSD

> Language: **English** · 中文總覽見 [../README.zh-TW.md](../README.zh-TW.md)

A fully scripted, EEGLAB-based resting-state EEG preprocessing pipeline, from EDF
to analysis-ready PSD (and optional band power). Function-based with centralized
config, deterministic seeds, version stamping, QC and logging.

## Requirements

- MATLAB R2020a or newer
- [EEGLAB](https://sccn.ucsd.edu/eeglab/) (2021+), installed and on the MATLAB path
- EEGLAB plugins: `clean_rawdata`, `ICLabel`, `firfilt`; `Picard` optional
  (falls back to `runica`)

> EEGLAB is **not** bundled. Run `addpath('/path/to/eeglab'); eeglab nogui;`
> first; the pipeline auto-detects EEGLAB's location.

## Two modes: automated vs. manual confirmation

First load EEGLAB in MATLAB:

```matlab
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
```

### Mode A — fully automated (`run_pipeline_auto`)

```matlab
run_pipeline_auto('/path/to/raw_edf');     % step1→step5, no human input
```

All ICA components are decided by ICLabel and epochs are rejected by rule; no
human intervention. Best for processing large batches quickly.

### Mode B — manual confirmation (`run_pipeline_manual`, **used for the manuscript**)

```matlab
run_pipeline_manual('/path/to/raw_edf');   % machine runs heavy steps; human confirms
```

The machine runs preprocessing and ICA, then **a human confirms two ambiguous
decisions** — the workflow used for the paper:

| Confirmation point | Script | What you confirm |
|--------------------|--------|------------------|
| ICA components | `review_ica_smart` | Only **mid-confidence** components (ICLabel max score inside the review window); high-confidence are auto-removed, low-confidence auto-kept. Shows topography / spectrum / activation; keys `r`=remove `k`=keep `s`=skip subject `q`=save & quit |
| Epochs | `review_epochs_smart` | Reviews epoch rejection for low-quality subjects |

> Both modes share identical step1–step5 logic and parameters; they differ only
> in whether ICA components and epochs are auto-confirmed or human-confirmed.
> Human decisions are written to `outputs/decisions/<ID>_ica_decisions.json` and
> reused on re-runs. Manual mode must be run in an **interactive** MATLAB session.

### Advanced: staged / custom runs

```matlab
% Manual mode with explicit start/stop steps
opts = struct('overwrite', false, 'skip_review', false, ...
              'start_from', 'step1', 'stop_after', 'step5');
master_pipeline([], opts);

% Re-run only some steps for one subject
master_pipeline({'sub01'}, struct('start_from','step3','stop_after','step4'));
```

## Steps

| Step | File | What it does |
|------|------|--------------|
| 1 | `step1_preprocess.m` | Load EDF, standardize channels (strip -LE, T3/4/5/6→T7/8/P7/P8), apply 10-20 locations, band-pass 0.5–70 Hz, multi-metric bad-channel detection (correlation / flatline / line-noise / kurtosis, ≥2 triggers), spherical interpolation |
| 2 | `step2_ica.m` | CleanLine 60 Hz (falls back to notch), rank-reduced ICA (Picard, falls back to runica), ICLabel classification |
| 2b | `step2b_apply_ica.m` | Apply removal decisions (Eye ≥ 0.70 / Muscle ≥ 0.80), average reference, ICA over-correction QC |
| 3 | `step3_epoch_reject.m` | 3 s epochs, three-layer rejection (autorej / amplitude 150 µV / jump 100 µV), sample 60 |
| 4 | `step4_export_psd.m` | Welch PSD (Hamming, no overlap), averaged across epochs → `psd_<ID>.mat` |
| 5 | `step5_export_bandpower.m` | (optional) absolute/relative band power from PSD → CSV |

## Input / output

- **Input**: `raw_edf/*.edf` (subject ID = EDF basename). Produce these from
  NeuroGuide with `../01_ng2edf`.
- **Output** (all under `outputs/`, gitignored):
  - `psd_mat/psd_<ID>.mat` — for `../03_features`
  - `band_power/` — optional MATLAB-side band power
  - `qc_pdf/`, `qc_png/`, `logs/`, `manifests/`, `methods/`

> Recommended: treat `../03_features` (Python) as the single source of truth for
> traditional band power, sharing the same PSD and band definition as FOOOF.
> `step5`'s MATLAB band power is optional / for cross-checking.

## Config

| File | Contents |
|------|----------|
| `config/config_paths.m` | Paths (CHAN_LOC points to bundled `chan_locs/`; EEGLAB_DIR auto-detected) |
| `config/config_thresholds.m` | All thresholds (filter, bad-channel, ICA, ICLabel, epoch, PSD, A/B/C/F grading) |
| `config/config_bands.m` | Band definitions |
| `config/config_pipeline.m` | Version, 19 channels, exclusion list, random seed |
| `config/config_user.m` | **Per-project override (gitignored)**. Copy `config_user.m.example` and edit |

## Other tools

- `master_qc_report.m` — 5-page QC PDF per subject (A/B/C/F grading)
- `master_methods_export.m` — auto-generate a COBIDAS-style Methods paragraph draft
- `population_qc_summary.m` — population-wide QC overview
- `install_picard.m` — install the Picard plugin (load EEGLAB first)
- `review_ica_smart.m` / `review_epochs_smart.m` — semi-automated manual review

## Reproducibility

- Deterministic per-subject seed (hash of subject ID + global seed)
- `version_stamp` records git hash + EEGLAB/plugin versions + timestamp
- Per-subject JSON manifest + per-step CSV log
- Note: `runica` is not fully deterministic; Picard is more stable but still
  depends on the BLAS implementation
