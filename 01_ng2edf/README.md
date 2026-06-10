# 01_ng2edf — NeuroGuide (.ng) → EDF

> Language: **English** · 中文總覽見 [../README.zh-TW.md](../README.zh-TW.md)

A reverse-engineered reader for the proprietary NeuroGuide `.ng` binary format,
with a standards-compliant EDF writer. NeuroGuide ships no public converter; this
script parses the header heuristically (length-prefixed ASCII strings +
fixed-offset integers) to recover the sample rate, sample count, channel count,
montage reference, and the per-channel float32 data block, keeping only the EEG
data channels.

## Requirements

- Python 3.8+
- `numpy` (see the top-level `requirements.txt`)

## Usage

```bash
python ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf

# Preview without writing
python ng2edf.py --input ./data/raw_ng --output ./data/raw_edf --dry-run

# Only process case IDs with a given prefix
python ng2edf.py -i ./data/raw_ng -o ./data/raw_edf --accept-prefix S
```

## Input layout

Each immediate sub-folder of `--input` is one subject and holds one `.ng`:

```
raw_ng/
├── S001=/        →  S001.edf
│   └── xxx.ng
├── S002=/        →  S002.edf
│   └── yyy.ng
└── ...
```

The case ID is the leading `[A-Za-z]+[0-9]+` token of the folder name, upper-cased
(`S001=` → `S001`, `A012X` → `A012`). If a folder has multiple `.ng` files, the
first (sorted) one is used.

## Channels

Default 19-channel 10-20 order:

```
FP1 FP2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 Fz Cz Pz
```

The montage reference (e.g. LinkEars) is appended as a suffix to channel names
(`FP1-LE`). The downstream MATLAB stage standardizes channels further (strips the
suffix, converts to modern 10-20 naming T7/T8/P7/P8).

## Privacy

The output EDF's *patient identification* field is set to the **de-identified
case ID** (from the folder name), **not** to any identifier that may be stored in
the NeuroGuide header. Ensure your input folder names are already de-identified.

## Notes

- Physical min/max use a single global range across all channels (one gain);
  keep this in mind for analyses sensitive to absolute amplitude.
- Parsing is heuristic reverse engineering — after conversion, spot-check a few
  files with an EDF viewer or MNE to confirm the waveforms look correct.
