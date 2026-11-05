# 01_ng2edf — NeuroGuide (.ng) → EDF

> Language: **English** · 中文總覽見 [../README.zh-TW.md](../README.zh-TW.md)

A reverse-engineered reader for the proprietary NeuroGuide `.ng` binary format,
with a standards-compliant EDF writer. NeuroGuide ships no public converter or
specification. The byte layout this reader follows is documented in
**[NG_FORMAT.md](NG_FORMAT.md)**.

The header is parsed structurally:

- The section offsets stored in the file locate the float32 data and the tail
  section.
- Section B gives the sample rate, the sample count and the electrode code of
  every data block.
- Only the EEG blocks (electrode codes 0..18) are kept.
- Consistency checks stop the conversion of a file whose layout does not match,
  instead of writing misaligned data.

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

# Also write <ID>_events.tsv with the event markers of each file
python ng2edf.py -i ./data/raw_ng -o ./data/raw_edf --events
```

| Option | Effect |
|---|---|
| `--events` | Also write `<ID>_events.tsv` next to each EDF (see [Event markers](#event-markers)). |
| `--lenient` | Convert files that fail the file-size or average-channel check, with a warning, instead of reporting them as errors. |

The exit code is 0 if every file converted and 2 if any file failed. Failed
files are listed as `ERROR` lines, and the batch continues.

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

EEG channels are selected by their NeuroGuide electrode code (0..18). They are
written in the 19-channel 10-20 order:

```
FP1 FP2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 Fz Cz Pz
```

The file's two pseudo-channels are not written to the EDF:

- the linked-ears reference (code 55), which is all zero;
- the average channel (code 56), which is used only as an integrity check.

The montage reference (e.g. LinkEars) is appended as a suffix to channel names
(`FP1-LE`). The downstream MATLAB stage standardizes channels further (strips the
suffix, converts to modern 10-20 naming T7/T8/P7/P8).

A channel that NeuroGuide stored as NaN (e.g. a removed electrode) is written
as 0, with a warning. Mark it as a bad channel downstream.

## Checks

Every file is checked against the layout in [NG_FORMAT.md](NG_FORMAT.md) §10.

- **Always an error**:
  - a missing magic;
  - inconsistent section offsets;
  - a section A or B that does not end at the next section offset;
  - a file truncated inside the EEG data.
- **Error, or a warning with `--lenient`**:
  - a file size other than tail offset + 5 × samples;
  - an average channel (code 56) that does not match the mean of the EEG
    channels. This would mean the channel table was misread.
- **Warning only**: unknown electrode codes, an unknown montage name,
  unexpected preamble constants, non-zero padding, and unrecognised tail
  records.

## Event markers

NeuroGuide stores event markers in the tail section of the file. `--events`
writes one row per non-empty tail record to `<ID>_events.tsv`. The file is
written even when a recording has no events, in which case it has a header
row only.

| Column | Meaning |
|---|---|
| `record` | index of the 6-byte tail record |
| `tail_offset`, `file_offset` | byte offset of the event code, from the start of the tail section and of the file |
| `code` | event code byte (only `2` has been seen) |
| `kind` | `event`, or `unknown` for a record with other non-zero bytes |
| `time_s_span`, `time_s_legacy`, `time_s_sample` | the event time in seconds under each candidate time base |
| `raw_hex` | the 6 record bytes |

**The time base of the tail records is not yet known.** The table therefore
reports the time under all three candidates. At 300 Hz they are record/250,
record/200 and record/300 s. Pick one only once it has been confirmed (see
[NG_FORMAT.md](NG_FORMAT.md) §8). `time_s_legacy` corresponds to the marker
positions of earlier converter versions, give or take about 20 ms. Event
markers are not written into the EDF.

## Privacy

- The output EDF's *patient identification* field is set to the
  **de-identified case ID** (from the folder name). It is **not** set to any
  identifier that may be stored in the NeuroGuide header. Ensure your input
  folder names are already de-identified.
- The NeuroGuide subject-information section can contain names, dates of birth
  and comments. The converter never decodes, prints or writes those strings.
- The converter reads the age field but never writes it. The EDF start
  date/time are taken from the NeuroGuide date/time of test.
- Treat the EDF files and the event tables as participant data. Both are
  excluded by the repository's `.gitignore` (`*.edf`, `*_events.tsv`).

## Tests

The tests use synthetic `.ng` files only (`tests/ngsynth.py` builds them in a
temporary folder). They include byte-for-byte checks of the EDF output against
the previous converter version. Run them from the repository root:

```bash
python -m unittest discover -s 01_ng2edf/tests -v
```

## Notes

- Physical min/max use a single global range across all channels (one gain);
  keep this in mind for analyses sensitive to absolute amplitude.
- The format is reverse engineered and has been verified on recordings from one
  acquisition setup (DSI-24, 19 channels, 300 Hz, linked ears). For other
  setups, check the warnings, and spot-check a few converted files with an EDF
  viewer or MNE.
