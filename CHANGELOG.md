# Changelog

All notable changes to this repository are listed here.

## [2.0.1] - 2026-11-21

### Changed
- `04_manuscript_sleepmed_r1/analysis/`: code tidied, with comments, docstrings and
  printed headings in English. Outputs are written to `outputs/analysis/`. The
  numerical results are unchanged.
- `04_manuscript_sleepmed_r1/README.md`: the preprocessing description now notes
  the visual check of flagged components and retained epochs, as in the manuscript.
- README (EN, zh-TW): shortened; the authorship section notes help from Codex
  with the Python refactoring.
- `.gitignore`: ignores `04_manuscript_sleepmed_r1/analysis/local_settings.json`.

## [2.0.0] - 2026-11-18

### Added
- `04_manuscript_sleepmed_r1/`: study-specific code for the revised *Sleep Medicine*
  manuscript — the preprocessing pipeline as run for the analysed data
  (`matlab_pipeline_v2.3`, including per-epoch and surface-Laplacian spectral
  export) and the Python code for the principal component analysis and
  clustering of the ISI items and its validation (`analysis/`).
  No participant data, identifiers or outputs are included.
- `04_manuscript_sleepmed_r1/README.zh-TW.md` and Traditional Chinese READMEs for
  each stage.

### Changed
- README (EN, zh-TW): new section "Authorship and de-identification";
  laboratory website updated to https://sleepaffect.org.
- `04_manuscript_sleepmed_r1/analysis/`: limited to the principal component
  analysis and clustering; statistical analyses and the code for tables and
  figures are not included.
- `01_ng2edf`: the NeuroGuide header is parsed structurally, following the new
  byte-level description in `01_ng2edf/NG_FORMAT.md`:
  - The data and tail offsets are read from the preamble instead of scanning
    for 0xFF padding.
  - The channel count and electrode codes are read from section B.
  - EEG channels are selected by electrode code (0..18) instead of by "first
    non-zero blocks".
  - The date and time of test come from their fixed slots in section A.
  - The EDF output is byte-identical to the previous version for files that
    converted before. The console output differs only by new warnings, such
    as the splice warning below.
- `01_ng2edf`: the header field earlier read as a "version" is the subject's
  age. It is now named `age` internally and is still not written to any output.

### Added
- `01_ng2edf`: integrity checks. The average channel (code 56) must equal the
  mean of the EEG channels (skipped when an EEG channel is entirely NaN), and
  the file size must equal tail offset + 5 × samples. A file that fails is
  reported and skipped. `--lenient` converts it with a warning.
- `01_ng2edf`: the tail section is decoded as one byte and one u32 event word
  per sample. The only event word seen, `0x00020000`, marks a splice: the
  recording is discontinuous after that sample, most likely because
  acquisition was paused. Splices are not user-placed Event Markers.
  - Every file with a splice gets a console warning with the splice times,
    because analysis windows must not span them. The EDF output is
    unchanged (plain EDF, no annotations).
  - `--events` writes `<ID>_events.tsv` with the columns `sample`, `time_s`,
    `word_hex`, `kind` (`splice` or `unknown`) and `note`.
  - Non-zero per-sample bytes in the tail, unrecognised event words and
    splices off the 6-sample block grid give warnings.
  - An earlier draft of this change read the tail as 6-byte records with
    three candidate time bases. That reading was wrong and has been removed,
    together with `event_time_s()` and `EVENT_TIME_BASES`. NG_FORMAT.md §11
    explains how the old readings relate to the correct one.
- `01_ng2edf/NG_FORMAT.md`: byte-level format description with field tables,
  invariants, open questions and privacy notes.
- `01_ng2edf/tests/`: unit tests on synthetic `.ng` files (stdlib `unittest`).
- `.gitignore`: `*_events.tsv`.

### Fixed
- `01_ng2edf`: a file whose first data byte is 0xFF is read from the correct
  offset.
- `01_ng2edf`: all-NaN channels are written as 0 with a warning instead of as
  garbage values.
