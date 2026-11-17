# Code for the Sleep Medicine manuscript (revision)

Study-specific code for the revised manuscript on insomnia symptom severity and
waking periodic beta activity (manuscript SLEEP-D-26-02563). The folders
`01_ng2edf`, `02_preprocess_eeglab` and `03_features` in this repository are the
general, reusable pipeline; this folder holds the code as it was run for the
manuscript.

**No participant data are included.** Raw recordings, questionnaire data,
participant identifiers and all analysis outputs stay with the research team, as
required by the ethics approval. The code therefore cannot be re-run without
access to the data, but every step reported in the manuscript can be read here.

```
04_manuscript_sleepmed_r1/
├── matlab_pipeline_v2.3/   # MATLAB + EEGLAB: preprocessing as run for the analysed data
└── analysis/               # Python: spectral parameterisation, statistics, tables, figures
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
| `step5_export_epoch_psd.m` | Power spectra for each epoch (epoch-level vigilance analysis) |
| `step6_export_psd_laplacian.m` | Spectra after a spherical-spline surface Laplacian |
| `validation_runica_vs_picard.m`, `validation_test_retest.m` | Checks of the ICA algorithm and of run-to-run reproducibility |
| `config/` | Paths, pipeline settings, thresholds, frequency bands |
| `CHANGELOG_v2.3.md` | Changes from the previous pipeline version |

`config/config_paths.m` is a template: data locations are read from the
environment variables `EEG_RAW_EDF_DIR`, `EEG_RAW_NG_DIR` and `EEGLAB_DIR`.
Participants excluded before preprocessing are listed in a local file,
`config/excluded_subjects.txt`, which is not part of the repository.

## `analysis/` — statistics, tables and figures

`run_r1_all.py` runs every step in order and records versions, input hashes and
the status of each step. Scripts can also be run one at a time.

| Script | Content | Manuscript |
|---|---|---|
| `task_b_fooof_full.py` | Spectral parameterisation (1–45 Hz) | Methods 2.4 |
| `task_f_recompute_relpower.py` | Conventional band power | Methods 2.6 |
| `task_e_bootstrap_stability.py`, `task_e2_bootstrap_selectK.py` | Bootstrap stability of the cluster solutions | Supplementary Table S1 |
| `task_r1w_pc_interpretation.py` | Principal components of the ISI items: parallel analysis, stability, external correlates, association with EEG measures | Results 3.1–3.2, Figure 1, Table 2, Table S1 |
| `task_r1b_k234_comparison.py` | Validity indices for two to six clusters; group contrasts, equivalence test, Bayes factors | Results 3.3, Table 2, Table S1 |
| `task_r1e_clustering_robustness.py` | Nine clustering methods and their agreement | Table S1 |
| `task_r1a_severity_form.py`, `task_r1g_followup.py`, `task_r1u_shape_robust.py` | Form of the association with severity: ISI bands, spline and step models, selection-corrected tests, robust standard errors, inclusion threshold | Results 3.3 and 3.6, Figures 2 and 6, Table S2 |
| `task_r1i_mood_adjustment.py` | Adjustment for depressive and anxiety symptoms | Results 3.3, Table 2 |
| `task_r1h_beta_characterisation.py`, `task_r1n_beta_subbands.py` | Frequency-wise comparison, scalp distribution, beta sub-bands | Results 3.4, Figures 3–4, Table S3 |
| `task_r1c_sensitivity_emg_vigilance.py`, `task_r1d_preproc_qc.py`, `task_r1t_surface_laplacian.py`, `task_r1s_vigilance_epochs.py` | Muscle activity, preprocessing quality, surface Laplacian, epoch-level vigilance | Results 3.5, Figure 5, Table S4 |
| `task_r1p_picard_comparison.py`, `task_r1j_manual_vs_auto.py`, `task_r1k_recruitment_flags.py` | Comparison of preprocessing runs; sample before and after the audit of recruitment records | Table S4 |
| `task_r1q_hrv_recompute.py`, `task_r1m_hrv_convergence.py` | Heart rate variability (exploratory) | Results 3.5, Table S5 |
| `task_r1f_pc_biomarkers.py`, `task_r1v_sleep_onset_check.py` | Additional analyses not reported in the manuscript | – |
| `make_r1_report_figures.py`, `make_r1_report_tables.py` | Intermediate report tables and figures | – |
| `sub_table1.py`, `sub_table2.py`, `sub_fig1.py` … `sub_fig6.py` | Tables and figures of the main text | Tables 1–2, Figures 1–6 |
| `sub_tableS1.py` … `sub_tableS5.py`, `sub_figS1.py`, `sub_figS2.py` | Supplementary tables and figures | Tables S1–S5, Figures S1–S2 |
| `_dynamic_n.py`, `_r1_common.py`, `_r1_spectra.py`, `_r1_submission.py`, `apa_docx_helper.py` | Shared code: sample definition, data loading, effect sizes, spectra, table and figure formatting | – |

### Inputs

- A SQLite database of questionnaire data (ISI items and total, PSQI, BDI-II,
  BAI, DBAS, age, sex, education); the columns read are listed in
  `_r1_common.py` and `_dynamic_n.py`.
- Spectra exported by `matlab_pipeline_v2.3` (`psd_<ID>.mat`, per-epoch spectra,
  Laplacian spectra).
- Electrocardiogram recordings for the heart rate variability analysis.

Locations and the identifiers of excluded participants are read from
`local_settings.json`, which is kept on the local machine;
`local_settings.example.json` shows the format. Environment variables
(`PSYML_DB`, `EEG_PIPELINE_ROOT`, …) take precedence.

### Environment

Python 3.12; package versions are in `analysis/requirements.txt`.
