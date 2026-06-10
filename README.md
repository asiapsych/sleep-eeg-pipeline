# Resting-State EEG Processing Pipeline (NeuroGuide → Features)

**English** | [中文](README.zh-TW.md)

A reproducible pipeline that takes raw NeuroGuide EEG recordings all the way to
analysis-ready **traditional band power** and **periodic / aperiodic (FOOOF)**
feature tables.

This repository contains **method code only** — no participant data. Downstream
subtyping, clinical measures and statistics are part of each study's private
analysis and are intentionally outside this open-source scope.

---

## Pipeline overview

A three-stage linear workflow, run in order:

```
  Raw .ng (NeuroGuide)
        │
        ▼
┌─────────────────────────┐
│ 01_ng2edf  (Python)     │   Reverse-engineered NeuroGuide reader → standard EDF
└─────────────────────────┘
        │  *.edf
        ▼
┌─────────────────────────┐
│ 02_preprocess_eeglab    │   filter → bad-channel detect/interp → CleanLine →
│   (MATLAB + EEGLAB)     │   ICA (ICLabel) → avg ref → epoch → reject → PSD export
└─────────────────────────┘
        │  psd_<ID>.mat
        ▼
┌─────────────────────────┐
│ 03_features  (Python)   │   Traditional band power (absolute / relative / ratios)
│                         │   + FOOOF aperiodic exponent/offset + periodic peaks
└─────────────────────────┘
        │
        ▼
   Feature tables (CSV)  ──►  (downstream statistics — not in this repo)
```

| Stage | Language | Input | Output |
|-------|----------|-------|--------|
| `01_ng2edf` | Python | folder of `.ng` | `.edf` |
| `02_preprocess_eeglab` | MATLAB + EEGLAB | `.edf` | `psd_<ID>.mat`, band-power CSV |
| `03_features` | Python | `psd_*.mat` | feature CSVs (traditional power + FOOOF) |

See each stage's `README.md` for full detail.

---

## Repository structure

```
.
├── README.md                 # this file (English)
├── README.zh-TW.md           # 中文版
├── LICENSE                   # BSD 3-Clause
├── requirements.txt          # pinned Python dependencies
├── .gitignore                # excludes all data + caches
├── 01_ng2edf/                # Python — NeuroGuide .ng → EDF
│   ├── ng2edf.py
│   └── README.md
├── 02_preprocess_eeglab/     # MATLAB — EDF → preprocessing → PSD (+ band power)
│   ├── run_pipeline_auto.m       # fully automated entry
│   ├── run_pipeline_manual.m     # human-in-the-loop entry (used for the manuscript)
│   ├── master_pipeline.m         # orchestrator (step1 … step5)
│   ├── step1_preprocess.m … step5_export_bandpower.m
│   ├── review_ica_smart.m / review_epochs_smart.m   # manual confirmation
│   ├── config/                   # paths, thresholds, bands, pipeline version
│   ├── helpers/                  # shared functions (bad-channel, rank, RNG, QC …)
│   ├── chan_locs/                # bundled 10-20 channel-location file
│   └── README.md
└── 03_features/              # Python — PSD → traditional power + FOOOF
    ├── compute_features.py
    └── README.md
```

---

## Requirements

### Python (stages 1 & 3)

```bash
pip install -r requirements.txt
```

- Stage 1 (conversion) needs only `numpy`.
- Stage 3 (features) needs `numpy / scipy / pandas / fooof`.
- FOOOF 1.x requires `numpy < 2` (pinned in `requirements.txt`).

### MATLAB (stage 2)

- MATLAB R2020a or newer.
- [EEGLAB](https://sccn.ucsd.edu/eeglab/) (2021 or newer), installed separately and added to the MATLAB path.
- EEGLAB plugins: `clean_rawdata`, `ICLabel`, `firfilt`; `Picard` optional (falls back to `runica`).
- EEGLAB is **not** bundled here. Add it yourself (`addpath('/path/to/eeglab'); eeglab nogui;`); the pipeline auto-detects its location.

---

## Quick start

```bash
# 1) NeuroGuide → EDF
python 01_ng2edf/ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf
```

```matlab
% 2) MATLAB preprocessing (inside an interactive MATLAB session)
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
run_pipeline_manual('/path/to/raw_edf')   % human-confirmed (manuscript mode)
% or run_pipeline_auto('/path/to/raw_edf') % fully automated
```

```bash
# 3) PSD → features
python 03_features/compute_features.py \
    --psd-dir /path/to/outputs/psd_mat --out-dir /path/to/features
```

---

## Preprocessing: two modes (automated vs. manual confirmation)

Stage 2 offers two equivalent one-command entry points that share identical
step1–step5 logic and parameters; they differ only in how the *ambiguous*
decisions are made:

- **`run_pipeline_auto`** — every ICA component and epoch is decided
  automatically (ICLabel + rule-based rejection). Best for fast batch processing.
- **`run_pipeline_manual`** *(used for the manuscript)* — the machine runs the
  heavy steps, then a human confirms two ambiguous decisions: mid-confidence
  ICA components (`review_ica_smart`, showing topography / spectrum / activation)
  and epoch rejection for low-quality subjects (`review_epochs_smart`).
  High/low-confidence cases are auto-handled; human decisions are saved to
  `outputs/decisions/<ID>_ica_decisions.json` and reused on re-runs.

---

## Design principles

- **Reproducible** — fixed random seeds; a deterministic per-subject seed
  (hash of subject ID + global seed); a version stamp (git hash + tool versions);
  full per-subject JSON manifests and per-step CSV logs.
- **Single source of truth for power** — traditional band power and FOOOF are
  both computed in stage 3 from the *same* PSD with the *same* band definition,
  so band edges never disagree between scripts.
- **Methods export** — `master_methods_export` (stage 2) auto-generates a
  COBIDAS-style Methods paragraph draft from the actual settings.

### Frequency bands (literature standard)

- **Primary bands** (continuous, non-overlapping; the denominator of relative
  power): Delta 1–4, Theta 4–8, Alpha 8–13, Beta 13–30, Gamma 30–45 Hz.
- **Derived sub-bands** (reported as a share of the same total, not in the
  denominator): Alpha1 8–10, Alpha2 10–13, Beta1 13–18, Beta2 18–24,
  Beta3 24–30, HighBeta 30–40 Hz.
- Relative power = band absolute power / sum of the 5 primary bands
  (= a single continuous 1–45 Hz integral) × 100.

---

## Data privacy

No raw EEG files, participant-level clinical data, database files, identifier
mappings, or study-specific subject codes are included in this repository.
Access to participant-level data is restricted by institutional ethics approval,
participant consent, and privacy regulations.

The pipeline is privacy-aware by construction: stage 1 writes a **de-identified
case code** (derived from the input folder name) into the EDF patient field and
never propagates any identifier that may be stored in the NeuroGuide header.

---

## License

This repository contains de-identified and generalized research code derived
from an internal EEG analysis pipeline developed for the associated manuscript.
The code is released under the **BSD 3-Clause License** (see `LICENSE`).

This public release was prepared by extracting, de-identifying, and generalizing
the reusable analysis components from the laboratory's internal study-specific
pipeline.

---

## Provenance and AI assistance

This public repository is a **de-identified and generalized rewrite** of the
laboratory's internal, study-specific analysis pipeline. The original code
contained participant-level file paths, identifier mapping tables, and database
couplings that cannot be shared publicly. The extraction, de-identification, and
generalization of the reusable components into this release were performed with
the assistance of **Claude (Claude Opus 4.8, Anthropic)**. All released code was
reviewed and approved by the authors, who are responsible for its correctness.

---

## Citation

The released code supports reproducibility of the EEG conversion, preprocessing,
PSD extraction, and feature-generation workflow. Please cite the associated
manuscript and this repository when reusing or adapting the code.

**Suggested Data/Code Availability statement (for the manuscript):**

> The de-identified and generalized analysis code used for EEG conversion,
> preprocessing, PSD extraction, and feature generation is available at
> https://github.com/YuCHii26/sleep-eeg-pipeline. Raw EEG recordings,
> participant-level clinical data, database files, and identifier mappings are
> not publicly available because the ethics approval and participant consent for
> this study did not include public sharing of participant-level data and
> because the data are subject to institutional privacy restrictions.

---

## Maintainers

**Sleep and Affective Science Lab** (睡眠與情緒科學研究室)
Division of Clinical Psychology, Department of Psychology, Asia University
Website: https://chang-lab-psychology-website-one.vercel.app/

## Acknowledgements

Built on [EEGLAB](https://sccn.ucsd.edu/eeglab/) (preprocessing/ICA), the
[ICLabel](https://github.com/sccn/ICLabel) and `clean_rawdata` plugins, and
[FOOOF / specparam](https://fooof-tools.github.io/fooof/) (spectral
parameterization).
