# NeuroGuide `.ng` file format (reverse-engineered)

This is a byte-level description of the NeuroGuide `.ng` EEG file, as read by
[`ng2edf.py`](ng2edf.py). NeuroGuide (Applied Neuroscience, Inc.) publishes no
specification, so everything here comes from reverse engineering.

**Evidence base.** The layout below was checked byte-for-byte on 30
NeuroGuide-produced recordings from a single acquisition setup (DSI-24 amplifier,
19 channels of the 10-20 system, 300 Hz, linked-ears montage). The header
constants and the tail section (§8) were also checked on the structural
metadata of 384 recordings from the same setup, and the meaning of the event
word was checked against the EEG of those 384 recordings. Statements marked
*invariant* held in every one of those files. They may not hold for other
amplifiers, montages, sample rates or NeuroGuide versions. The converter
checks the ones it depends on.

**Confidence.** H = structure and meaning established. M = structure
established, meaning likely. L = structure established, meaning unknown.

Contents: [1 Conventions](#1-conventions) ·
[2 File layout](#2-file-layout) · [3 Preamble](#3-preamble-33-bytes) ·
[4 Section A](#4-section-a-subject-information) ·
[5 Section B](#5-section-b-signal-and-montage-description) ·
[6 Data](#6-data-section) · [7 Gap](#7-gap-21-bytes) ·
[8 Tail and events](#8-tail-section-per-sample-bytes-and-event-words) ·
[9 Invariants](#9-invariants) · [10 What ng2edf checks](#10-what-ng2edf-checks) ·
[11 Earlier readings](#11-earlier-readings-of-the-format) ·
[12 Open questions](#12-open-questions) · [13 Privacy](#13-privacy) ·
[14 References](#14-references)

---

## 1. Conventions

- All integers are little-endian. `u8`/`u16`/`u32` are unsigned and `i16` is
  signed. `f32` is IEEE-754 single precision.
- **string**: a `u32` length *L* that includes the terminating NUL, followed by
  *L* bytes. NeuroGuide writes an **empty string as L = 1 plus a single NUL**,
  so an empty string is the 5 bytes `01 00 00 00 00`. The converter also
  accepts L = 0 with no bytes. Text is 8-bit; Windows-1252 is assumed, and only
  ASCII has been seen.
- Offsets written `0x..` are absolute. "B+k" means k bytes after the start of
  section B.

## 2. File layout

```
0          33              off_B             off_data                       off_tail                     EOF
| preamble | A: subject info | B: signal/montage | f32 data [n_ch][n] | 21 B | u8 [n] | u32 words [n] |
```

`n` is the number of samples per channel and `n_ch` the number of data blocks
(both from section B). The preamble stores the four section offsets. A reader
must take them from there. It must not assume a fixed header length, because
the header grows with every non-empty subject string, and it must not scan for
0xFF bytes.

Header length, for reference:
`off_data = 33 + len(A) + len(B)`, where

- `len(A) = Σᵢ (4 + Lᵢ) + 1`, summed over the 13 strings, plus the single u8
  field; and
- `len(B) = 8 + 2·n_ch + 3 + Σ_tables (2 + 4·rows) + (4 + L_name) + 2 + 2·n₃`.

For example, the synthetic test file has a 569-byte header (empty strings
except a 5-character age, a 10-character date, a 3-character label and an
8-character time, with 21 channels, two 19-row montage tables, the name
"LinkEars" and a 110-entry table 3).

## 3. Preamble (33 bytes)

| Offset | Type | Value seen | Meaning | Conf. |
|---|---|---|---|---|
| 0x00 | char[10] | `NeuroGuide` | magic, no length prefix, no NUL | H |
| 0x0a | u8 | 0 | reserved / NUL after magic | L |
| 0x0b | u16 × 3 | 2, 4, 8 | unknown constants (format id? sizes of u16/u32/f64?) | L |
| 0x11 | u32 | 33 | offset of section A | H |
| 0x15 | u32 | | offset of section B (= end of section A) | H |
| 0x19 | u32 | | **offset of the data section** (= end of section B = header length) | H |
| 0x1d | u32 | | **offset of the tail section** = off_data + 4·n_ch·n + 21 | H |

## 4. Section A: subject information

Section A holds 13 strings, with a single `u8` field between strings #11 and
#12:

```
str#1 … str#11, u8, str#12, str#13
```

| # | Role | Format | Conf. | Notes |
|---|---|---|---|---|
| 1, 2 | unknown | free text | – | empty in the files examined |
| 3 | **age** in years | `%.2f`, e.g. `34.00` | H | Earlier readers called this a "version". It is not: NeuroGuide stores the age in each `.ng` for its age-regressed norms. Age is a quasi-identifier. |
| 4–6 | unknown | free text | – | empty in the files examined |
| 7 | **date of test** | `MM/DD/YYYY` | H | |
| 8–11 | unknown | free text | – | empty in the files examined |
| (u8) | unknown | 1 in all files examined | L | Possibly a radio-button setting (gender, handedness, eyes closed/open?) |
| 12 | **subject label** | free text | H (position), L (which field) | NeuroGuide's Name, Subject ID or EEG ID field. It can hold a real name. |
| 13 | **time of test** | `HH:MM:SS` | H | Whether this is the start or the save time is unknown. |

The NeuroGuide *Subject Information* window has 13 fields: Name, Subject ID,
Date of Birth, Age, Gender, Handedness, EEG ID, Date of Test, Time of Test,
Technician, Clinician, Medication and Comments. The stored order is not the
window order. Only slots 3, 7, 12 and 13 have been identified, because the other
fields were empty in every file examined. **Any of the unidentified slots may
hold a name, a date of birth or free-text clinical comments.**

How the strings are delimited. An earlier reading took the first byte of
every empty string (`01`) to be a per-string "flag" (u8 = 1) followed by a
length-0 string, giving the pattern `F F U F F F U F F F F F U`. The flagged
positions are exactly the empty strings (plus #12, which follows the u8
field). Under that reading a non-empty string would also need a flag, but
NeuroGuide writes none. Both readings give the same strings for empty
fields. The reading adopted here (uniform strings, empty = length 1) also
covers non-empty strings in any slot. `ng2edf.py` still accepts the older
form: before each string it allows an optional `0x01` byte, and the u8 field
may be missing. Each slot is resolved structurally. A layout is accepted
only if all 13 strings are well formed and the last one ends exactly at
`off_B`.

## 5. Section B: signal and montage description

| Position | Type | Value seen | Meaning | Conf. |
|---|---|---|---|---|
| B+0 | u16 | 300 | sample rate, Hz | H |
| B+2 | u32 | | n, samples per channel. It sets the data and tail sizes. | H |
| B+6 | u16 | 21 | n_ch, number of f32 data blocks | H |
| B+8 | u16 × n_ch | 0..18, 55, 56 | **electrode code** of each data block, in data order (table below) | H (0..18), M-H (55, 56) |
| next | u16 | 55 | reference electrode code (linked ears) | M-H |
| next | u8 | 0x47 | unknown | L |
| next | u16 + (u16 active, u16 ref) × rows | 19 × (i, 55) | montage table 1: channel *i* referenced to code 55 | H (structure) |
| next | u16 + (u16, u16) × rows | same as table 1 | montage table 2 | M (structure), L (purpose) |
| next | string | `LinkEars` | montage name | H |
| next | u16 n₃ + i16 × n₃ | 110 × −1 | "table 3", unknown. These 220 `0xFF` bytes are what older readers scanned for. | L |
| end | | | equals `off_data` | H |

Electrode codes:

| Code | Channel | Code | Channel | Code | Channel |
|---|---|---|---|---|---|
| 0 | FP1 | 7 | P4 | 14 | T5 |
| 1 | FP2 | 8 | O1 | 15 | T6 |
| 2 | F3 | 9 | O2 | 16 | Fz |
| 3 | F4 | 10 | F7 | 17 | Cz |
| 4 | C3 | 11 | F8 | 18 | Pz |
| 5 | C4 | 12 | T3 | **55** | linked-ears reference (LE). Its block is all zero (LE re LE). |
| 6 | P3 | 13 | T4 | **56** | average of the EEG channels (§6) |

Codes above 56 have not been seen. NeuroGuide's montage editor supports more
scalp and auxiliary channels, so other recordings may use more codes.

## 6. Data section

At `off_data`: `n_ch` blocks of `n` × f32 each, channel-major (all samples of
block 0, then block 1, …), in the order of the electrode-code table.

- EEG blocks (codes 0..18) are in µV against the montage reference. µV is
  inferred from the amplitudes; no unit field has been found.
- The block with **code 55** is all `0.0`.
- The block with **code 56** equals `float32( (Σ_c x_c as float64) / 19 )`,
  summed over the 19 EEG blocks. The match is bit-exact. It is NaN wherever an
  EEG channel is NaN.
- **NaN channels**: a channel can be stored as NaN at every sample, with the bit
  pattern `0xFFC00000` (the x86 default quiet NaN). No partial-NaN channels and
  no ±inf have been seen.
- The first data byte can be `0xFF`, since it is the low byte of a float. A
  reader that skips "all 0xFF bytes" before the data starts one byte late.

## 7. Gap (21 bytes)

`off_tail − (off_data + 4·n_ch·n)` = **21** bytes, all zero in every file
examined. 21 equals n_ch, so this might be one byte per channel (for example a
bad-channel flag), but that is a guess.

## 8. Tail section: per-sample bytes and event words

From `off_tail` to EOF: **5·n bytes**, one byte and one `u32` word per sample.

| Tail offset | Type | Content seen | Meaning | Conf. |
|---|---|---|---|---|
| 0 … n−1 | u8 × n | all zero | one byte per sample; meaning unknown (an edit-selection or rejection state?) | H (structure), L (meaning) |
| n … 5n−1 | u32 × n | zero except at events | **event word of sample *i***, at tail offset n + 4*i* | H |

**Events.** An event is a non-zero event word. Its index is the sample index
*s* (0-based), so its time is **s / fs** seconds from the first sample. The
only value seen is **`0x00020000`** (bytes `00 00 02 00`). Its single
non-zero byte is at tail offset p = n + 4s + 2, so s = (p − n − 2) / 4.

Every `0x00020000` word seen is on the **last sample of a 6-sample block**
(s ≡ 5 mod 6, a 20 ms block at 300 Hz), and n has always been a multiple of 6.
NeuroGuide apparently appends the amplifier data in 6-sample blocks.

**`0x00020000` marks a splice.** The recording is discontinuous between
sample *s* and sample *s* + 1. The evidence, from 384 recordings:

- At 86 of the 87 words found, all channels jump at the same time, exactly
  between *s* and *s* + 1. The jump is in the top 1 % of the file's
  sample-to-sample differences, and usually the largest jump in the file.
  This held for all 21 words in the 30 byte-checked recordings.
- No comparable jump was found at a block boundary without a word.
- A recording that contains a splice is shorter than the same participant's
  other recordings of the same protocol (median difference about 9 s). This
  fits a paused interval that is missing from the file.
- After high-pass filtering, a brief broadband transient follows the splice
  and lasts about 1 s.

The most likely cause is an acquisition **Pause → Continue**. (NeuroGuide's
version history introduces Pause/Continue together with event markers.) An
amplifier reset or an impedance check during acquisition could leave the same
trace, so this is not settled (§12). The word is **not** a user-placed Event
Marker: it has a single value, and its position is quantised to the
acquisition block. How user-placed Event Markers are stored is unknown; none
are present in the files examined.

**For analysis**, treat every splice as a segment boundary:

- The EEG before and after a splice is not contiguous. If the splice is a
  pause, the file time after it also lags the wall-clock time by the paused
  interval. As far as is known, the file does not store that interval.
- Do not let epochs, FFT/Welch windows, filters or spindle, slow-wave and
  artifact detectors span a splice. Process the two sides separately, and
  discard at least 1 s on either side of the splice; this also covers the
  transient.
- The EDF written by `ng2edf.py` is a plain, continuous EDF and does not mark
  splices. Take them from the converter's warning or the `--events` table.
  The discontinuity lies between `time_s` and `time_s + 1/fs`.

**What `ng2edf.py` does.** It decodes both arrays. Every file with a splice
gets one warning that lists the splice times. `--events` writes one row per
non-zero word, with the columns `sample`, `time_s`, `word_hex`, `kind`
(`splice` for `0x00020000`, `unknown` for any other value) and `note`.
Non-zero per-sample bytes, unknown words and splices off the 6-sample grid
also give warnings. The per-sample bytes do not change the EDF.

An earlier version of this description read the tail as 6-byte records and
listed three candidate time bases. That reading was wrong (§11).

## 9. Invariants

These held in every file examined. The converter checks those marked ✔
(§10).

1. ✔ Magic `NeuroGuide`. ✔ The offsets increase: 33 ≤ off_A < off_B < off_data < off_tail ≤ file size.
2. off_A = 33.
3. ✔ Section A ends exactly at off_B and section B ends exactly at off_data.
4. ✔ off_tail ≥ off_data + 4·n_ch·n. The difference was 21 in every file (✔ warning only).
5. ✔ File size = off_tail + 5·n (n per-sample bytes, then n u32 event words).
6. n is a multiple of 6.
7. Constants: preamble `2, 4, 8` (✔ warning only), byte 0x0a = 0, u8 field = 1,
   u8 after the reference code = 0x47, montage tables 19 × (i, 55) twice,
   montage name `LinkEars`, table 3 = 110 × −1.
8. Channel table `0..18, 55, 56`. The block with code 55 is all zero. ✔ The
   block with code 56 is the mean of the EEG blocks.
9. Non-finite values occur only as whole channels of `0xFFC00000`.
10. The gap and the per-sample byte array are all zero (✔ warning only).
    Every non-zero event word is `0x00020000`, at a sample s ≡ 5 (mod 6) (✔
    other values or positions give a warning).

## 10. What ng2edf checks

| Check | On failure |
|---|---|
| magic, offset order, offsets within the file | error |
| section A: 13 well-formed strings ending at off_B | error (a fit that deviates from the known layout gives a warning) |
| section B: complete, ending at off_data | error |
| float data fits before off_tail | error |
| file size = off_tail + 5·n | error, or a warning with `--lenient` |
| code-56 block = mean of the EEG blocks (within 1e-3 µV + 1e-5 relative), on samples where all EEG is finite. Skipped when an EEG channel is entirely non-finite. | error, or a warning with `--lenient`. This catches a misread channel table. |
| at least one EEG code (0..18), no duplicates | error |
| unknown electrode codes | warning; the block is ignored |
| unknown montage name | warning; labels get no reference suffix |
| date / time of test not in the expected format | warning; the EDF falls back to `01.01.00` / `00.00.00` |
| preamble constants, gap size/content | warning |
| tail: non-zero per-sample bytes (with their count), unrecognised event words, a splice off the 6-sample grid, non-zero bytes in an incomplete last word | warning |
| splices (event word `0x00020000`) | one warning per file, with the splice times |
| non-finite EEG samples | warning; the samples are set to 0 (mark the channel as bad downstream) |

EEG channels are selected **by electrode code (0..18)**, not by position or by
testing for non-zero data. They are labelled `FP1 … Pz` plus a
reference suffix (`-LE` for `LinkEars`).

## 11. Earlier readings of the format

These are listed so that results from older readers can be related to this
description.

- **"220 bytes of 0xFF padding"**: table 3 (u16 count 110, then 110 × −1). Its
  length depends on the count, and the data offset is stored at 0x19.
- **"u32 n_map = 21, then codes 1..18, 55, 56, 55"**: the same bytes read 2
  bytes off. The correct reading is u16 n_ch = 21, then u16 codes 0..18, 55,
  56, then the u16 reference code 55.
- **"`47 13 00 00 00` then pairs `(i, 0x37)`"**: the u8 0x47, then the u16
  row count 19 of montage table 1, then its rows `(i, 55)`.
- **"22 data blocks, block 21 = markers"** (the earlier converter): there are
  21 f32 blocks. The "block 21" reader took one u32 per sample, starting 21
  bytes before the tail section, i.e. n + 21 bytes before the event words.
  - Its sample index is s_old = s + ⌊(n + 23)/4⌋. That is about a quarter of
    the recording too late.
  - It could not see events with s ≥ (3n − 23)/4, i.e. in roughly the last
    quarter of a file.
  - Its two "codes", 512 and 2²⁵, are the word `0x00020000` read through a
    misaligned u32 frame.
  - To convert its marker positions: s = s_old − ⌊(n + 23)/4⌋.
- **"6-byte tail records, event code 2 in byte 4"** (an earlier version of
  this description):
  - The 6-byte period is an artefact of the block grid. Splices sit at
    s = 6m + 5, so their `02` byte is at tail offset n + 24m + 22, which is
    ≡ 4 (mod 6) because n is a multiple of 6.
  - "Code 2" is byte 2 of the word `0x00020000`.
  - The three candidate time bases derived from it (record/250, record/200
    and record/300 s at 300 Hz) are all wrong.
- **"NG version" (a string such as `34.00`)**: this is the age field (§4).
- **"F-flag before strings"**: the length byte of an empty string (§4).

## 12. Open questions

1. **Preamble constants** `2, 4, 8` at 0x0b and the byte at 0x0a.
2. **The u8 field before string #12** (= 1). This is what earlier readings
   saw as "flag bytes". Which setting is it?
3. **Roles of strings 1, 2, 4–6 and 8–11**, and whether #12 is the Name,
   Subject ID or EEG ID field.
4. **Time of test**: is it the recording start or the save time?
5. **`0x47` after the reference code**: montage type, flags or a setting?
6. **Why two montage tables?** Perhaps one for recording and one for display.
   They have always been identical.
7. **Table 3** (110 × −1): per-electrode-code values that were never set? It
   could hold impedances, bad-channel marks or a code-to-position map.
8. **The 21 zero bytes** between the data and the tail: one byte per data
   block?
9. **The per-sample byte array** (§8). NeuroGuide says `.ng` stores "edited
   EEG selections", so it may hold a selection or rejection state per
   sample. It was all zero in the files examined.
10. **Which operator action writes `0x00020000`.** Pause → Continue is the
    most likely; an amplifier reset or an impedance check are alternatives.
    What its two u16 halves mean (a type and a flag?) is a guess.
11. **Where user-placed Event Markers are stored** (Eyes Open, Unspecified …).
    None occur in the files examined: no other event-word value was seen.
    They may use other values in the same word array, or be stored elsewhere.
12. **Where the eyes-open / eyes-closed condition is stored.** Candidates are
    the u8 field before #12, `0x47`, or the tail.
13. **Other setups**: other amplifiers, sample rates, channel counts, montages
    (average, Laplacian) and NeuroGuide versions.

Most of these could be settled by saving test recordings from NeuroGuide
with **dummy data**. For example: press Pause and Continue at known times
during a recording, fill each Subject Information field with a distinct
dummy text, toggle eyes open/closed, place event markers at known times,
make an edit selection or mark a bad channel, then compare the saved files
byte by byte. Use only dummy data; never real participant data.

## 13. Privacy

- **Section A can contain direct identifiers**: a name, a date of birth, IDs,
  technician or clinician names, medication and free-text comments. It also
  holds the **age** and the **date and time of test**, which are
  quasi-identifiers.
- `ng2edf.py` **never copies section A into its outputs**. It locates the
  strings by their lengths only, and it decodes only the age, date and time
  slots after checking their fixed formats on the raw bytes. No other string
  is decoded, printed, logged or written. Error and warning messages contain
  only offsets, lengths, counts, sample indices and times from the start of
  the recording.
- The **age is not written** to the EDF or anywhere else.
- The **EDF patient field** is the de-identified case ID taken from the input
  folder name.
- The **EDF start date and time** are copied from the date and time of test,
  as the EDF header requires them. The `--events` table contains sample
  indices, times relative to the start of the recording and event-word
  values, but no header content.
- Do **not** publish or share `.ng` files, and do not publish hex dumps of
  their headers. Treat converted EDF files and event tables as participant
  data too.

## 14. References

- Applied Neuroscience, Inc. *NeuroGuide Manual*,
  <https://www.appliedneuroscience.com/PDFs/NeuroGuide_Manual.pdf>: the
  Subject Information window; "Save the Edited EEG Selections and patient
  information in NeuroGuide Format (*.ng)"; the montage stored in the `.ng`
  file; the Event Marker chapter.
- Applied Neuroscience, Inc. *History of NeuroGuide*,
  <https://www.appliedneuroscience.com/PDFs/History_of_NG.pdf>: version
  history, event markers and DSI-24 support.
- No public `.ng` specification or other open-source reader was found.
