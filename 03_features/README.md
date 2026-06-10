# 03_features — PSD → traditional power + FOOOF (periodic / aperiodic)

> Language: **English** · 中文總覽見 [../README.zh-TW.md](../README.zh-TW.md)

From the `psd_<ID>.mat` files exported by stage 02, compute two feature families
from the **same PSD with the same band definition**:

1. **Traditional band power** — absolute / relative power + common ratios (TBR, ABR, …)
2. **Periodic / aperiodic decomposition (FOOOF / specparam)** — aperiodic exponent
   and offset, plus per-band periodic peaks (center frequency / power / bandwidth)

This script is **feature extraction only**. Group assignment, clinical measures
and statistical tests belong to each study's private analysis and are not here.

## Requirements

- Python 3.8+
- `numpy / scipy / pandas / fooof` (see the top-level `requirements.txt`;
  FOOOF 1.x requires numpy < 2)

## Usage

```bash
python compute_features.py --psd-dir /path/to/outputs/psd_mat --out-dir /path/to/features

# Skip specific subjects; customize the fit frequency range
python compute_features.py --psd-dir ./psd_mat --out-dir ./features \
    --exclude sub01 sub02 --freq-min 1 --freq-max 45
```

## Input

`psd_<ID>.mat` (stage 02 `step4_export_psd` output), each containing:

| Variable | Shape | Meaning |
|----------|-------|---------|
| `freqs` | `[n_freq]` | frequency axis (Hz) |
| `psd` | `[n_ch × n_freq]` | power spectral density |
| `ch_names` | `[n_ch]` | channel labels |
| `n_ep` | scalar (optional) | number of epochs averaged |

## Output (CSV)

| File | Contents |
|------|----------|
| `traditional_power_subject.csv` | per subject × absolute/relative power per band (whole-brain + 5 regions) + ratios |
| `specparam_channel.csv` | per subject × per channel: exponent / offset / peaks / R² |
| `specparam_subject.csv` | whole-brain / regional averages of the specparam features |
| `fit_quality.csv` | FOOOF fit R² / error per subject |

## Frequency bands (literature standard)

- **Primary bands** (continuous, non-overlapping; denominator of relative power):
  Delta 1–4, Theta 4–8, Alpha 8–13, Beta 13–30, Gamma 30–45 Hz
- **Derived sub-bands** (reported as a share of the total, not in the denominator):
  Alpha1 8–10, Alpha2 10–13, Beta1 13–18, Beta2 18–24, Beta3 24–30, HighBeta 30–40 Hz
- Relative power = band absolute power / sum of the 5 primary bands
  (= a single continuous 1–45 Hz integral) × 100
- Ratios (TBR, ABR, …) use absolute power (the total cancels out)

## Regions (19 channels, modern 10-20)

| Region | Channels |
|--------|----------|
| FRONTAL | FP1 FP2 F3 F4 F7 F8 FZ |
| CENTRAL | C3 C4 CZ |
| PARIETAL | P3 P4 PZ |
| TEMPORAL | T7 T8 P7 P8 |
| OCCIPITAL | O1 O2 |

> If the input PSD uses legacy temporal naming (T3/T4/T5/T6), it is normalized to
> modern naming (T7/T8/P7/P8) on load, so region computation is unaffected.

## FOOOF settings

`freq_range=(1,45)`, `aperiodic_mode="fixed"`, `peak_width_limits=(1,12)`,
`max_n_peaks=8`, `min_peak_height=0.1`, `peak_threshold=2.0`.
