# Changelog

All notable changes to this repository are listed here.

## [Unreleased]

### Changed
- `01_ng2edf`: the NeuroGuide header is parsed structurally, following the new
  byte-level description in `01_ng2edf/NG_FORMAT.md`:
  - The data and tail offsets are read from the preamble instead of scanning
    for 0xFF padding.
  - The channel count and electrode codes are read from section B.
  - EEG channels are selected by electrode code (0..18) instead of by "first
    non-zero blocks".
  - The date and time of test come from their fixed slots in section A.
  - The EDF output is byte-identical to the previous version for files that
    converted before.
- `01_ng2edf`: the header field earlier read as a "version" is the subject's
  age. It is now named `age` internally and is still not written to any output.

### Added
- `01_ng2edf`: integrity checks. The average channel (code 56) must equal the
  mean of the EEG channels (skipped when an EEG channel is entirely NaN), and
  the file size must equal tail offset + 5 × samples. A file that fails is
  reported and skipped. `--lenient` converts it with a warning.
- `01_ng2edf`: `--events` writes `<ID>_events.tsv` with the event markers of
  the tail section. Each row has the record index, byte offsets, code and the
  time under each of three candidate time bases; the time base is not yet
  known.
- `01_ng2edf/NG_FORMAT.md`: byte-level format description with field tables,
  invariants, open questions and privacy notes.
- `01_ng2edf/tests/`: unit tests on synthetic `.ng` files (stdlib `unittest`).
- `.gitignore`: `*_events.tsv`.

### Fixed
- `01_ng2edf`: a file whose first data byte is 0xFF is read from the correct
  offset.
- `01_ng2edf`: all-NaN channels are written as 0 with a warning instead of as
  garbage values.
