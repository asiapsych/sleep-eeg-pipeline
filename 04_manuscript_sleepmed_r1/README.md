# Code for the Sleep Medicine manuscript (revision)

> Language: **English** · [繁體中文](README.zh-TW.md)

Study-specific code for the revised manuscript on insomnia symptom severity and
waking periodic beta activity (manuscript SLEEP-D-26-02563). The folders
`01_ng2edf`, `02_preprocess_eeglab` and `03_features` in this repository are the
general, reusable pipeline; this folder holds the preprocessing as it was run
for the manuscript and the code for the principal component analysis and
clustering of the Insomnia Severity Index (ISI) items.

**No participant data are included.** Raw recordings, questionnaire data,
participant identifiers and all analysis outputs stay with the research team, as
required by the ethics approval. The code therefore cannot be re-run without
access to the data, but each step can be read here.

```
04_manuscript_sleepmed_r1/
├── matlab_pipeline_v2.3/   # MATLAB + EEGLAB: preprocessing as run for the analysed data
└── analysis/               # Python: principal component analysis and clustering of the ISI items
```

Comments and docstrings in these two folders are in Traditional Chinese.

## `matlab_pipeline_v2.3/` — preprocessing

Fully automated, with one set of fixed thresholds for all participants
(`config/config_thresholds.m`): band-pass filter, bad-channel detection and
interpolation, CleanLine, ICA with PICARD, ICLabel-based component removal,
average reference, epoching and epoch rejection, then spectral export.

| File | Purpose |
|---|---|
| `run_all_v23.m` | Entry point: steps 1–5 for all recordings |
| `run_step6_v23.m` | Entry point: step 6 only |
| `step1_preprocess.m` … `step3_epoch_reject.m` | Filtering, channels, ICA, epochs |
| `step4_export_psd.m` | Power spectra averaged over 60 epochs |
| `step5_export_epoch_psd.m` | Power spectra for each epoch |
| `step6_export_psd_laplacian.m` | Spectra after a spherical-spline surface Laplacian |
| `validation_runica_vs_picard.m`, `validation_test_retest.m` | Checks of the ICA algorithm and of run-to-run reproducibility |
| `config/` | Paths, pipeline settings, thresholds, frequency bands |
| `CHANGELOG_v2.3.md` | Changes from the previous pipeline version |

`config/config_paths.m` is a template: data locations are read from the
environment variables `EEG_RAW_EDF_DIR`, `EEG_RAW_NG_DIR` and `EEGLAB_DIR`.
Participants excluded before preprocessing are listed in a local file,
`config/excluded_subjects.txt`, which is not part of the repository.

## `analysis/` — principal component analysis and clustering

All procedures use the seven ISI items of participants with an ISI total of 9 or
above. Scripts are run one at a time.

| Script | Content |
|---|---|
| `task_r1w_pc_interpretation.py` | Principal component analysis of the seven ISI items: loadings, explained variance, parallel analysis (2,000 permuted data sets) and bootstrap congruence of the loadings (2,000 resamples) |
| `task_e_bootstrap_stability.py` | Loader of the clustering sample; bootstrap silhouette and predict-back adjusted Rand index for the three-cluster solution (reference for `task_e2`) |
| `task_e2_bootstrap_selectK.py` | Predict-back bootstrap stability (1,000 resamples) for two to six clusters |
| `task_r1b_k234_comparison.py` | Internal validity indices (silhouette, Calinski–Harabasz, Davies–Bouldin) for two to six clusters, with the bootstrap stability from `task_e2` |
| `task_r1e_clustering_robustness.py` | Nine clustering methods (K-means, four Gaussian mixture models, three hierarchical methods, latent class analysis) for two and three clusters, their agreement, and the gap statistic |
| `_dynamic_n.py`, `_r1_common.py` | Shared code: sample definition and data loading |

Run `task_e2_bootstrap_selectK.py` before `task_r1b_k234_comparison.py` and
`task_r1e_clustering_robustness.py`, which read its bootstrap results.

### Inputs

- A SQLite database of questionnaire data (ISI items and total, age, sex,
  education); the columns read are listed in `_r1_common.py` and `_dynamic_n.py`.
- Spectral feature tables (periodic and aperiodic parameters and relative band
  power; see stage 03), which define the participants with usable EEG.

The database location and the identifiers of excluded participants are read
from `local_settings.json`, which is kept on the local machine;
`local_settings.example.json` shows the format. The environment variable
`PSYML_DB` takes precedence.

### Environment

Python 3.12; package versions are in `analysis/requirements.txt`.
