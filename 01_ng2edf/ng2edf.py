#!/usr/bin/env python3
"""
ng2edf.py — NeuroGuide (.ng) → EDF converter
=============================================

Reverse-engineered reader for the proprietary NeuroGuide ``.ng`` binary
format, with a standards-compliant EDF writer. Scans a folder of subject
sub-folders (one ``.ng`` each) and writes one EDF per subject.

The file is parsed structurally, following the byte-level description in
``NG_FORMAT.md`` (same folder): the preamble stores the offsets of the
float32 data and of the tail section, section B gives the sample rate, the
sample count and the electrode code of every data block, and only the EEG
blocks (electrode codes 0..18, FP1 … Pz) are written to the EDF. Consistency
checks (average channel, file size) stop the conversion of a file whose
layout does not match, instead of writing misaligned data.

Privacy note
------------
The output EDF's *patient identification* field is set to the
**de-identified case ID** derived from the input folder name — NOT to any
identifier that may be stored inside the NeuroGuide header. Make sure your
input folder names are already de-identified (e.g. ``S001=`` → ``S001``).

The NeuroGuide subject-information strings (which can hold a name, a date
of birth, comments …) are located but never decoded, printed or written.
Only the date and time of test are copied, into the EDF start date/time
fields. The age field is read (``info['age']``) but never written anywhere.

Usage
-----
    python ng2edf.py --input  /path/to/raw_ng \
                     --output /path/to/raw_edf

    python ng2edf.py --input ./data/raw_ng --output ./data/raw_edf --dry-run

    # also write <ID>_events.tsv with the event markers of each file
    python ng2edf.py --input ./data/raw_ng --output ./data/raw_edf --events

Folder naming
-------------
Each immediate sub-folder of ``--input`` is treated as one subject. The
case ID is the leading ``[A-Za-z]+[0-9]+`` token of the folder name,
upper-cased (``s003=`` → ``S003``, ``A012X`` → ``A012``).
Use ``--accept-prefix`` to restrict which prefixes are processed.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import struct
import sys

import numpy as np

# Standard 19-channel 10-20 order as stored by NeuroGuide (electrode codes
# 0..18 in the channel table of section B).
NG_CHANNEL_NAMES_19 = [
    'FP1', 'FP2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2',
    'F7', 'F8', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz',
]


# ── NeuroGuide format constants (see NG_FORMAT.md) ───────────────────────

NG_MAGIC = b'NeuroGuide'
PREAMBLE_SIZE = 33                   # magic, u8, 3 x u16, 4 x u32 offsets
PREAMBLE_CONSTS = (2, 4, 8)          # u16 x 3 at 0x0b, meaning unknown

# Section A: 13 length-prefixed strings, with one u8 field between strings
# #11 and #12 (see _parse_section_a).
SECTION_A_STRINGS = 13
SECTION_A_BYTE_BEFORE = 12           # the u8 field precedes string #12
LEGACY_FLAG = 0x01                   # older reading: optional u8 flag per string
SLOT_AGE, SLOT_DATE, SLOT_TIME = 3, 7, 13

N_EEG_CODES = 19                     # electrode codes 0..18 = FP1 … Pz
CODE_LINKED_EARS = 55                # reference pseudo-channel (all zero)
CODE_AVERAGE = 56                    # average of the EEG channels

GAP_BYTES = 21                       # between float data and tail
TAIL_BYTES_PER_SAMPLE = 5            # tail section = 5 x n_samples bytes
TAIL_RECORD_SIZE = 6                 # the tail is a sequence of 6-byte records
TAIL_EVENT_BYTE = 4                  # event code byte inside a record

# The time base of the tail records is NOT known. Candidate record rates,
# as (records, samples): 'span' = the records cover the whole recording
# (5 per 6 samples, 250/s at 300 Hz); 'legacy' = the reading of the earlier
# converter (2 per 3 samples, 200/s at 300 Hz); 'sample' = one record per
# sample (300/s at 300 Hz). The events table reports all three.
EVENT_TIME_BASES = {'span': (5, 6), 'legacy': (2, 3), 'sample': (1, 1)}

# Montage names that become a channel-label suffix ('LinkEars' → '-LE').
KNOWN_REFERENCES = ('LinkEars', 'LinkedEars', 'Avg', 'Average',
                    'AvgRef', 'CzRef', 'A1Ref', 'A2Ref')

# Format checks for the machine-written strings (applied to raw bytes, so
# nothing else is ever decoded).
_RE_AGE = re.compile(rb'\d{1,3}\.\d{2}')
_RE_DATE = re.compile(rb'\d{2}/\d{2}/\d{4}')
_RE_TIME = re.compile(rb'\d{2}:\d{2}:\d{2}')


class NGFormatError(ValueError):
    """The bytes do not match the known .ng layout (see NG_FORMAT.md)."""


# ── NeuroGuide header parsing ────────────────────────────────────────────

def _unpack(fmt, data, pos, end, what):
    """struct.unpack_from that refuses to read past ``end``."""
    size = struct.calcsize(fmt)
    if pos + size > end:
        raise NGFormatError(f'{what} at byte {pos} runs past the end of its section')
    return struct.unpack_from(fmt, data, pos)


def _string_length(data, pos, end):
    """Declared length of the string at ``pos``, or None if it cannot be one.

    A string is a u32 length (including the trailing NUL) followed by that
    many bytes. NeuroGuide writes an empty string as length 1 (a lone NUL);
    length 0 (no bytes) is accepted as well. Only the length and the
    terminating NUL are inspected, never the text.
    """
    if pos + 4 > end:
        return None
    length = struct.unpack_from('<I', data, pos)[0]
    if pos + 4 + length > end:
        return None
    if length and data[pos + 4 + length - 1] != 0:
        return None
    return length


def _parse_section_a(data, start, end):
    """Locate the 13 subject-information strings in ``data[start:end]``.

    Returns ``(slots, deviates, ambiguous)``. ``slots`` lists
    ``(offset, length)`` for each string (offset of its u32 length), or is
    None when no layout fits.

    Known layout: strings #1-#11, one u8 field, strings #12-#13. NeuroGuide
    writes an empty string as length 1 plus a NUL, so an empty string is
    the 5 bytes ``01 00 00 00 00``. Earlier readings took the first of these
    bytes for a per-string flag (followed by a length-0 string). That reading
    gives the same strings for empty fields, but it also expects a flag
    before non-empty strings, where NeuroGuide writes none.

    Every slot is resolved structurally, so that a file whose bytes deviate
    from the known layout is still read correctly when possible: before
    each string an optional 0x01 byte (the older flag reading) is allowed,
    and the u8 field before #12 may be missing. A layout is accepted only if
    all strings are well formed and the last one ends exactly at the
    section-B offset. The known layout is tried first; ``deviates`` is True
    if it did not fit, ``ambiguous`` if two layouts with different strings
    both fit. Layouts that differ only in how an empty string is spelled
    count as the same.
    """
    layouts = {}

    def walk(slot, pos, acc, deviations):
        if len(layouts) > 1:
            return
        if slot == SECTION_A_STRINGS:
            if pos == end:
                key = tuple((o + 4 + n, max(n - 1, 0)) for o, n in acc)
                layouts.setdefault(key, (list(acc), deviations))
            return
        if slot + 1 == SECTION_A_BYTE_BEFORE:
            options = [(pos + 1, 0), (pos, 1)] if pos < end else [(pos, 1)]
        elif pos < end and data[pos] == LEGACY_FLAG:
            options = [(pos, 0), (pos + 1, 1)]
        else:
            options = [(pos, 0)]
        for p, extra in options:
            length = _string_length(data, p, end)
            if length is None:
                continue
            acc.append((p, length))
            walk(slot + 1, p + 4 + length, acc, deviations + extra)
            acc.pop()

    walk(0, start, [], 0)
    if not layouts:
        return None, False, False
    slots, deviations = next(iter(layouts.values()))
    return slots, deviations > 0, len(layouts) > 1


def _formatted_field(data, slot, pattern):
    """Value of a machine-written string ('' if empty, None if malformed)."""
    offset, length = slot
    if length == 0:
        return ''
    raw = data[offset + 4:offset + 4 + length - 1]
    if pattern.fullmatch(raw):
        return raw.decode('ascii')
    return None


def parse_header(data):
    """Parse the preamble, section A (string positions only) and section B.

    ``data`` is the complete file content. Returns a dict; raises
    :class:`NGFormatError` when the header does not match the layout in
    ``NG_FORMAT.md``. Free-text subject strings are never decoded.
    """
    file_size = len(data)
    if data[:len(NG_MAGIC)] != NG_MAGIC:
        raise NGFormatError('missing "NeuroGuide" magic')
    if file_size < PREAMBLE_SIZE:
        raise NGFormatError('file is shorter than the preamble')

    byte_0a = data[10]
    consts = struct.unpack_from('<3H', data, 0x0b)
    off_a, off_b, off_data, off_tail = struct.unpack_from('<4I', data, 0x11)
    if not (PREAMBLE_SIZE <= off_a < off_b < off_data < off_tail <= file_size):
        raise NGFormatError(
            f'section offsets out of order or beyond the end of the file '
            f'(A={off_a}, B={off_b}, data={off_data}, tail={off_tail}, '
            f'file size={file_size})')

    # ── section A: subject information (13 strings) ──────────────────
    slots, deviates, ambiguous = _parse_section_a(data, off_a, off_b)
    if slots is None:
        raise NGFormatError(
            f'section A (bytes {off_a}..{off_b}) does not hold the expected '
            f'13 length-prefixed strings')
    date = _formatted_field(data, slots[SLOT_DATE - 1], _RE_DATE)
    time_str = _formatted_field(data, slots[SLOT_TIME - 1], _RE_TIME)
    age = _formatted_field(data, slots[SLOT_AGE - 1], _RE_AGE)

    # ── section B: signal and montage description ────────────────────
    pos = off_b
    sample_rate, total_samples, n_ch = _unpack('<HIH', data, pos, off_data,
                                               'signal description')
    pos += 8
    if sample_rate == 0 or total_samples == 0:
        raise NGFormatError(
            f"Invalid sample_rate={sample_rate} or total_samples={total_samples}")
    if n_ch == 0:
        raise NGFormatError('channel table is empty')
    codes = list(_unpack(f'<{n_ch}H', data, pos, off_data, 'channel table'))
    pos += 2 * n_ch
    ref_code, byte_after_ref = _unpack('<HB', data, pos, off_data, 'reference code')
    pos += 3
    for table in (1, 2):
        rows, = _unpack('<H', data, pos, off_data, f'montage table {table}')
        pos += 2 + 4 * rows
    name_len = _string_length(data, pos, off_data)
    if name_len is None:
        raise NGFormatError(f'montage name at byte {pos} is not a valid string')
    montage = data[pos + 4:pos + 4 + max(name_len - 1, 0)]
    pos += 4 + name_len
    t3_count, = _unpack('<H', data, pos, off_data, 'table 3')
    pos += 2 + 2 * t3_count
    if pos != off_data:
        raise NGFormatError(
            f'section B ends at byte {pos}, but the data offset is {off_data}')

    known = {name.encode('ascii'): name for name in KNOWN_REFERENCES}
    return {
        'file_size': file_size,
        'byte_0a': byte_0a,
        'consts': consts,
        'subject_offset': off_a,
        'signal_offset': off_b,
        'data_offset': off_data,
        'tail_offset': off_tail,
        'section_a_deviates': deviates,
        'section_a_ambiguous': ambiguous,
        'date': date,
        'time': time_str,
        'age': age,
        'sample_rate': sample_rate,
        'total_samples': total_samples,
        'n_channels': n_ch,
        'channel_codes': codes,
        'ref_code': ref_code,
        'byte_after_ref': byte_after_ref,
        'reference': known.get(montage, ''),
    }


# ── Tail section: event markers ──────────────────────────────────────────

def parse_events(tail, tail_offset=0):
    """Return the non-empty 6-byte records of a tail section.

    ``tail`` is the tail section (bytes from the tail offset to the end of
    the file); ``tail_offset`` is its absolute position, used for
    ``file_offset``. A record whose only non-zero byte is byte 4 is an event
    (``kind='event'``, ``code`` = that byte); any other non-zero record is
    returned with ``kind='unknown'``. No time is attached: see
    :func:`event_time_s`.
    """
    n_rec = len(tail) // TAIL_RECORD_SIZE
    recs = np.frombuffer(tail[:n_rec * TAIL_RECORD_SIZE], dtype=np.uint8)
    recs = recs.reshape(n_rec, TAIL_RECORD_SIZE)
    events = []
    for r in np.flatnonzero(recs.any(axis=1)):
        rec = recs[r]
        others = np.delete(rec, TAIL_EVENT_BYTE)
        offset = int(r) * TAIL_RECORD_SIZE + TAIL_EVENT_BYTE
        events.append({
            'record': int(r),
            'tail_offset': offset,
            'file_offset': tail_offset + offset,
            'code': int(rec[TAIL_EVENT_BYTE]),
            'kind': 'unknown' if others.any() or not rec[TAIL_EVENT_BYTE] else 'event',
            'raw_hex': rec.tobytes().hex(),
        })
    return events


def event_time_s(event, sample_rate, time_base):
    """Time of ``event`` in seconds under an explicitly chosen time base.

    ``time_base`` is a key of ``EVENT_TIME_BASES`` ('span', 'legacy' or
    'sample'). The record rate of the tail section is unresolved (see
    NG_FORMAT.md), so there is deliberately no default.
    """
    if time_base not in EVENT_TIME_BASES:
        raise ValueError(f'time_base must be one of {sorted(EVENT_TIME_BASES)}')
    records, samples = EVENT_TIME_BASES[time_base]
    return event['record'] * samples / (records * sample_rate)


EVENT_TSV_COLUMNS = (['record', 'tail_offset', 'file_offset', 'code', 'kind']
                     + [f'time_s_{base}' for base in EVENT_TIME_BASES]
                     + ['raw_hex'])


def write_events_tsv(filepath, events, sample_rate):
    """Write the tail records as a tab-separated table (header line first)."""
    with open(filepath, 'w', encoding='ascii', newline='\n') as f:
        f.write('\t'.join(EVENT_TSV_COLUMNS) + '\n')
        for ev in events:
            row = [str(ev['record']), str(ev['tail_offset']),
                   str(ev['file_offset']), str(ev['code']), ev['kind']]
            row += [f'{event_time_s(ev, sample_rate, base):.6f}'
                    for base in EVENT_TIME_BASES]
            row.append(ev['raw_hex'])
            f.write('\t'.join(row) + '\n')


# ── Whole-file reader ────────────────────────────────────────────────────

def parse_ng(filepath, lenient=False):
    """Parse a NeuroGuide ``.ng`` file into a dict of header + channel data.

    Raises :class:`NGFormatError` (a ``ValueError``) when the file does not
    match the known layout. With ``lenient=True`` a wrong file size or an
    average channel that does not match the EEG only produce a warning.
    """
    with open(filepath, 'rb') as f:
        data = f.read()

    if data[0:10] != NG_MAGIC:
        raise NGFormatError(f"Not a NeuroGuide file: {filepath}")

    hdr = parse_header(data)
    warnings = []
    checks = {}

    def problem(msg):
        if not lenient:
            raise NGFormatError(f'{msg} (re-run with --lenient to convert anyway)')
        warnings.append(msg)

    if hdr['byte_0a'] != 0 or hdr['consts'] != PREAMBLE_CONSTS:
        warnings.append('unrecognised preamble constants; the layout may '
                        'differ from NG_FORMAT.md')
    if hdr['section_a_deviates']:
        warnings.append('section A deviates from the known layout; its strings '
                        'were located by search')
    if hdr['section_a_ambiguous']:
        warnings.append('section A has more than one consistent string layout; '
                        'date and time of test may be wrong')
    if hdr['date'] is None:
        warnings.append('date of test is not MM/DD/YYYY (EDF start date left at default)')
    if hdr['time'] is None:
        warnings.append('time of test is not HH:MM:SS (EDF start time left at default)')

    n = hdr['total_samples']
    codes = hdr['channel_codes']
    n_ch = len(codes)
    data_offset = hdr['data_offset']
    tail_offset = hdr['tail_offset']
    data_end = data_offset + 4 * n_ch * n
    if data_end > tail_offset:
        raise NGFormatError(
            f'the float data ({n_ch} x {n} samples) would end at byte {data_end}, '
            f'past the tail offset {tail_offset}')

    expected_size = tail_offset + TAIL_BYTES_PER_SAMPLE * n
    if len(data) != expected_size:
        problem(f'file size is {len(data)} bytes, expected {expected_size} '
                f'(tail offset + 5 x n_samples)')
    gap = data[data_end:tail_offset]
    if len(gap) != GAP_BYTES or any(gap):
        warnings.append(f'{len(gap)} bytes between data and tail section '
                        f'({GAP_BYTES} zero bytes expected)')

    # ── channel selection by electrode code ──────────────────────────
    eeg_idx = [i for i, c in enumerate(codes) if c < N_EEG_CODES]
    if not eeg_idx:
        raise NGFormatError(f"No EEG data channels found in: {filepath}")
    eeg_codes = [codes[i] for i in eeg_idx]
    if len(set(eeg_codes)) != len(eeg_codes):
        raise NGFormatError('the channel table lists an EEG electrode twice')
    unknown = sorted({c for c in codes
                      if c >= N_EEG_CODES and c not in (CODE_LINKED_EARS, CODE_AVERAGE)})
    if unknown:
        warnings.append(f'ignoring channel(s) with unknown electrode codes {unknown}')

    blocks = np.frombuffer(data[data_offset:data_end], dtype='<f4').reshape(n_ch, n)
    channels = [blocks[i] for i in eeg_idx]

    # ── integrity: the code-56 block is the mean of the EEG blocks ───
    if CODE_AVERAGE in codes:
        finite = np.ones(n, dtype=bool)
        total = np.zeros(n, dtype=np.float64)
        for ch in channels:
            finite &= np.isfinite(ch)
            total += ch
        if finite.any():
            avg = blocks[codes.index(CODE_AVERAGE)][finite].astype(np.float64)
            mean = total[finite] / len(channels)
            ok = np.abs(avg - mean) <= 1e-3 + 1e-5 * np.abs(mean)
            checks['average'] = 'ok' if ok.all() else 'mismatch'
            if not ok.all():
                problem('the average channel (code 56) does not match the mean '
                        'of the EEG channels; the channel table may be misread')
        else:
            checks['average'] = 'skipped (non-finite EEG channel)'
    else:
        checks['average'] = 'absent'

    if hdr['reference']:
        reference = hdr['reference']
    else:
        reference = ''
        warnings.append('montage name not recognised; channel labels get no '
                        'reference suffix')
    ref_suffix = ''
    if reference == 'LinkEars':
        ref_suffix = '-LE'
    elif reference:
        ref_suffix = f'-{reference}'

    labels = [NG_CHANNEL_NAMES_19[c] + ref_suffix for c in eeg_codes]

    # A channel stored as NaN (e.g. an electrode removed in NeuroGuide) cannot
    # be written to EDF. Set it to 0 and warn, so it is not written as garbage.
    for c in range(len(channels)):
        bad = ~np.isfinite(channels[c])
        if bad.any():
            channels[c] = np.where(bad, 0.0, channels[c]).astype(channels[c].dtype)
            warnings.append(f"{labels[c]} has {int(bad.sum())} non-finite samples "
                            f"(set to 0; mark it as a bad channel)")

    # ── tail section: event markers ──────────────────────────────────
    tail = data[tail_offset:tail_offset + TAIL_BYTES_PER_SAMPLE * n]
    events = parse_events(tail, tail_offset)
    n_unknown = sum(ev['kind'] != 'event' for ev in events)
    if n_unknown:
        warnings.append(f'{n_unknown} tail record(s) with unrecognised content '
                        f'(listed with --events)')
    if any(tail[len(tail) - len(tail) % TAIL_RECORD_SIZE:]):
        warnings.append('non-zero bytes after the last whole tail record')

    for msg in warnings:
        print(f"  WARNING: {msg}: {filepath}")

    return {
        'sample_rate': hdr['sample_rate'],
        'total_samples': n,
        'n_eeg': len(channels),
        'channels': channels,
        'date': hdr['date'] or '',
        'time': hdr['time'] or '',
        'age': hdr['age'] or '',          # never written to the EDF
        'reference': reference,
        'labels': labels,
        'channel_codes': eeg_codes,
        'data_offset': data_offset,
        'tail_offset': tail_offset,
        'events': events,
        'checks': checks,
        'warnings': warnings,
    }


# ── EDF writer ───────────────────────────────────────────────────────────

def write_edf(filepath, info, patient_id):
    """Write a 16-bit EDF. ``patient_id`` is the de-identified case ID."""
    channels = info['channels']
    labels = info['labels']
    sample_rate = info['sample_rate']
    total_samples = info['total_samples']
    n_channels = len(channels)
    samples_per_record = sample_rate
    n_records = math.ceil(total_samples / samples_per_record)

    global_min = min(float(ch.min()) for ch in channels)
    global_max = max(float(ch.max()) for ch in channels)
    phys_min = math.floor(global_min) - 1
    phys_max = math.ceil(global_max) + 1

    dig_min = -32768
    dig_max = 32767
    gain = (phys_max - phys_min) / (dig_max - dig_min)

    try:
        parts = info['date'].split('/')
        edf_date = f"{parts[1]}.{parts[0]}.{parts[2][-2:]}"
    except (IndexError, ValueError):
        edf_date = '01.01.00'

    edf_time = info['time'].replace(':', '.') if info['time'] else '00.00.00'

    header_bytes = 256 + n_channels * 256

    with open(filepath, 'wb') as f:
        f.write(b'0       ')
        # Patient identification = de-identified case ID (not the NG header).
        f.write(patient_id.ljust(80)[:80].encode('ascii', errors='replace'))
        f.write(b''.ljust(80))                                   # recording id
        f.write(edf_date.ljust(8)[:8].encode('ascii'))
        f.write(edf_time.ljust(8)[:8].encode('ascii'))
        f.write(str(header_bytes).ljust(8).encode('ascii'))
        f.write(b''.ljust(44))
        f.write(str(n_records).ljust(8).encode('ascii'))
        f.write(b'1'.ljust(8))                                   # record duration (s)
        f.write(str(n_channels).ljust(4).encode('ascii'))

        for label in labels:
            f.write(label.ljust(16)[:16].encode('ascii'))
        for _ in range(n_channels):
            f.write(b''.ljust(80))                               # transducer
        for _ in range(n_channels):
            f.write(b'uV'.ljust(8))                              # physical dim
        for _ in range(n_channels):
            f.write(str(phys_min).ljust(8)[:8].encode('ascii'))
        for _ in range(n_channels):
            f.write(str(phys_max).ljust(8)[:8].encode('ascii'))
        for _ in range(n_channels):
            f.write(str(dig_min).ljust(8).encode('ascii'))
        for _ in range(n_channels):
            f.write(str(dig_max).ljust(8).encode('ascii'))
        for _ in range(n_channels):
            f.write(b''.ljust(80))                               # prefiltering
        for _ in range(n_channels):
            f.write(str(samples_per_record).ljust(8).encode('ascii'))
        for _ in range(n_channels):
            f.write(b''.ljust(32))                               # reserved

        zero_pad = np.full(samples_per_record,
                           round((0.0 - phys_min) / gain + dig_min), dtype=np.int16)
        for rec in range(n_records):
            start = rec * samples_per_record
            end = min(start + samples_per_record, total_samples)
            actual = end - start
            for ch in range(n_channels):
                phys = channels[ch][start:end].astype(np.float64)
                digital = np.clip(
                    np.round((phys - phys_min) / gain + dig_min),
                    dig_min, dig_max,
                ).astype(np.int16)
                if actual < samples_per_record:
                    digital = np.concatenate([digital, zero_pad[actual:]])
                f.write(digital.tobytes())


# ── Folder scanning / case IDs ───────────────────────────────────────────

def extract_case_id(dirname):
    """Leading ``[A-Za-z]+[0-9]+`` token of a folder name, upper-cased."""
    m = re.match(r'^([A-Za-z]+\d+)', dirname)
    return m.group(1).upper() if m else None


def find_ng_cases(raw_ng_dir, accept_prefixes):
    """Return [(case_id, ng_path)] — first ``.ng`` (sorted) per sub-folder."""
    cases = []
    for entry in sorted(os.listdir(raw_ng_dir)):
        full = os.path.join(raw_ng_dir, entry)
        if not os.path.isdir(full):
            continue
        case_id = extract_case_id(entry)
        if case_id is None:
            continue
        if accept_prefixes and not any(case_id.startswith(p) for p in accept_prefixes):
            continue
        ng_files = sorted(f for f in os.listdir(full) if f.lower().endswith('.ng'))
        if ng_files:
            cases.append((case_id, os.path.join(full, ng_files[0])))
    return cases


# ── CLI ──────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Batch-convert NeuroGuide .ng files to EDF.')
    parser.add_argument('--input', '-i', required=True,
                        help='Folder of subject sub-folders, each with one .ng')
    parser.add_argument('--output', '-o', required=True,
                        help='Output folder for .edf files')
    parser.add_argument('--accept-prefix', action='append', default=None,
                        metavar='PREFIX',
                        help='Only process case IDs starting with PREFIX '
                             '(repeatable). Default: accept all.')
    parser.add_argument('--dry-run', action='store_true',
                        help='List planned conversions without writing files')
    parser.add_argument('--events', action='store_true',
                        help='Also write <ID>_events.tsv with the event markers '
                             'of each file (record index, byte offsets and the '
                             'candidate times; see NG_FORMAT.md)')
    parser.add_argument('--lenient', action='store_true',
                        help='Convert files that fail the file-size or '
                             'average-channel check (with a warning) instead '
                             'of skipping them')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.input):
        parser.error(f"input folder not found: {args.input}")
    os.makedirs(args.output, exist_ok=True)

    cases = find_ng_cases(args.input, args.accept_prefix)
    if not cases:
        print('No .ng files found.')
        return 1

    print(f"Found {len(cases)} subject(s).\n")
    if args.dry_run:
        for case_id, ng_path in cases:
            print(f"  {case_id}.edf  <-  {os.path.relpath(ng_path, args.input)}")
        print("\n--dry-run: nothing written.")
        return 0

    converted = failed = 0
    for case_id, ng_path in cases:
        edf_path = os.path.join(args.output, f'{case_id}.edf')
        try:
            info = parse_ng(ng_path, lenient=args.lenient)
            duration = info['total_samples'] / info['sample_rate']
            print(f"  {case_id}: {info['n_eeg']}ch | {info['sample_rate']}Hz | "
                  f"{info['total_samples']} samples ({duration:.0f}s) | "
                  f"ref={info['reference'] or 'n/a'}")
            write_edf(edf_path, info, patient_id=case_id)
            print(f"    -> {case_id}.edf")
            if args.events:
                tsv_name = f'{case_id}_events.tsv'
                write_events_tsv(os.path.join(args.output, tsv_name),
                                 info['events'], info['sample_rate'])
                print(f"    -> {tsv_name} ({len(info['events'])} record(s))")
            converted += 1
        except Exception as e:                                    # noqa: BLE001
            print(f"  ERROR {case_id} ({ng_path}): {e}")
            failed += 1

    print(f"\nDone: {converted} converted, {failed} failed.")
    return 0 if failed == 0 else 2


if __name__ == '__main__':
    sys.exit(main())
