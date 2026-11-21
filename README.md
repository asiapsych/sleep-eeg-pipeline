# Resting EEG processing and analysis

**English** | [中文](README.zh-TW.md)

Stages 1–3 convert NeuroGuide EEG recordings to EDF, preprocess the recordings
and extract spectral features. The repository contains code only. The *Sleep
Medicine* manuscript's preprocessing, principal component analysis and
clustering of the ISI items are in
[`04_manuscript_sleepmed_r1/`](04_manuscript_sleepmed_r1/README.md).

---

## Processing stages

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
├── CHANGELOG.md
├── 01_ng2edf/                # Python, NeuroGuide .ng to EDF
│   ├── ng2edf.py
│   ├── NG_FORMAT.md              # byte-level description of the .ng format
│   ├── tests/                    # unit tests on synthetic .ng files
│   └── README.md
├── 02_preprocess_eeglab/     # MATLAB, preprocessing and PSD export
│   ├── run_pipeline_auto.m       # fully automated entry
│   ├── run_pipeline_manual.m     # human-in-the-loop entry (used for the manuscript)
│   ├── master_pipeline.m         # orchestrator (step1 … step5)
│   ├── step1_preprocess.m … step5_export_bandpower.m
│   ├── review_ica_smart.m / review_epochs_smart.m   # manual confirmation
│   ├── config/                   # paths, thresholds, bands, pipeline version
│   ├── helpers/                  # shared functions (bad-channel, rank, RNG, QC …)
│   ├── chan_locs/                # bundled 10-20 channel-location file
│   └── README.md
├── 03_features/              # Python, band power and FOOOF from PSD
│   ├── compute_features.py
│   └── README.md
└── 04_manuscript_sleepmed_r1/   # study-specific code for the Sleep Medicine manuscript
    ├── matlab_pipeline_v2.3/     # preprocessing as run for the analysed data
    ├── analysis/                 # PCA and clustering of the ISI items (Python)
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
# 1) Convert NeuroGuide files to EDF   (add --events for a per-file table of recording splices)
python 01_ng2edf/ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf
```

Stage 1 prints a warning for every recording that contains a *splice*, a
discontinuity most likely left by an acquisition pause. The EDF does not mark
it, so keep epochs and filters from spanning it (see
[`01_ng2edf/README.md`](01_ng2edf/README.md#events-recording-splices)).

```matlab
% 2) MATLAB preprocessing (inside an interactive MATLAB session)
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
run_pipeline_manual('/path/to/raw_edf')   % human-confirmed (manuscript mode)
% or run_pipeline_auto('/path/to/raw_edf') % fully automated
```

```bash
# 3) Extract features from PSD files
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

- Each run uses fixed seeds, records the Git revision and tool versions, and
  writes JSON manifests and CSV logs. The seed for each recording is derived
  from its subject ID and the global seed.
- Stage 3 computes traditional band power and FOOOF from the same PSD using the
  same band definitions.
- `master_methods_export` (stage 2) auto-generates a
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

Stage 1 writes a de-identified case code from the input folder name into the EDF
patient field. It does not copy identifiers from the NeuroGuide header
(see [`01_ng2edf/NG_FORMAT.md`](01_ng2edf/NG_FORMAT.md#13-privacy)).

---

## License

This repository contains de-identified and generalized research code derived
from an internal EEG analysis pipeline developed for the associated manuscript.
The code is released under the **BSD 3-Clause License** (see `LICENSE`).

This public release was prepared by extracting, de-identifying, and generalizing
the reusable analysis components from the laboratory's internal study-specific
pipeline.

---

## Authorship, de-identification and AI assistance

The authors developed the research pipeline. Claude (Anthropic) assisted with
separating reusable code from the internal scripts and removing participant
identifiers, participant-level filenames and local paths before release.
Claude Opus 4.8 assisted with stages 01–03, and Claude Opus 5.5 assisted with
`04_manuscript_sleepmed_r1/`. Identifiers and paths are read from private local
settings or replaced with examples. The authors reviewed the de-identified code
before its initial release.

Codex (OpenAI) assisted with refactoring the Python PCA and clustering scripts.
We reran the scripts and compared their numerical outputs with the previous
version.

---

## Citation

The released code supports reproducibility of the EEG conversion, preprocessing,
PSD extraction, and feature-generation workflow. Please cite the associated
manuscript and this repository when reusing or adapting the code.

Archived version used for the revised manuscript: v2.0.2,
https://doi.org/10.5281/zenodo.24415381.

**Suggested Data/Code Availability statement (for the manuscript):**

> The de-identified and generalized analysis code used for EEG conversion,
> preprocessing, PSD extraction, and feature generation is available at
> https://github.com/sleepaffect/insomnia-resting-eeg. Raw EEG recordings,
> participant-level clinical data, database files, and identifier mappings are
> not publicly available because the ethics approval and participant consent for
> this study did not include public sharing of participant-level data and
> because the data are subject to institutional privacy restrictions.

---

## Maintainers

**Sleep and Affective Science Lab** (睡眠與情緒科學研究室)
Division of Clinical Psychology, Department of Psychology, Asia University
Website: https://sleepaffect.org

## Acknowledgements

Built on [EEGLAB](https://sccn.ucsd.edu/eeglab/) (preprocessing/ICA), the
[ICLabel](https://github.com/sccn/ICLabel) and `clean_rawdata` plugins, and
[FOOOF / specparam](https://fooof-tools.github.io/fooof/) (spectral
parameterization).
