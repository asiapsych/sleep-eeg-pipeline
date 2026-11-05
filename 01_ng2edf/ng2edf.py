#!/usr/bin/env python3
"""
ng2edf.py — NeuroGuide (.ng) → EDF converter
=============================================

Reverse-engineered reader for the proprietary NeuroGuide ``.ng`` binary
format, with a standards-compliant EDF writer. Scans a folder of subject
sub-folders (one ``.ng`` each) and writes one EDF per subject.

The NeuroGuide header is parsed heuristically (length-prefixed ASCII
strings + fixed-offset integers) to recover sample rate, sample count,
channel count, montage reference and the per-channel float32 data block.
Only the EEG data channels are kept.

Privacy note
------------
The output EDF's *patient identification* field is set to the
**de-identified case ID** derived from the input folder name — NOT to any
identifier that may be stored inside the NeuroGuide header. Make sure your
input folder names are already de-identified (e.g. ``S001=`` → ``S001``).

Usage
-----
    python ng2edf.py --input  /path/to/raw_ng \
                     --output /path/to/raw_edf

    python ng2edf.py --input ./data/raw_ng --output ./data/raw_edf --dry-run

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

# Standard 19-channel 10-20 order as stored by NeuroGuide.
NG_CHANNEL_NAMES_19 = [
    'FP1', 'FP2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2',
    'F7', 'F8', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz',
]


# ── NeuroGuide header parsing ────────────────────────────────────────────

def _scan_length_prefixed_strings(data, start, end):
    """Yield (offset, length, text) for length-prefixed ASCII strings."""
    strings = []
    pos = start
    while pos < end:
        if pos + 4 > len(data):
            break
        slen = struct.unpack_from('<I', data, pos)[0]
        if 2 <= slen <= 200 and pos + 4 + slen <= len(data):
            candidate = data[pos + 4:pos + 4 + slen]
            if all(32 <= b < 127 or b == 0 for b in candidate):
                text = candidate.rstrip(b'\x00').decode('ascii', errors='replace')
                if len(text) >= 1:
                    strings.append((pos, slen, text))
                    pos = pos + 4 + slen
                    continue
        pos += 1
    return strings


def parse_ng(filepath):
    """Parse a NeuroGuide ``.ng`` file into a dict of header + channel data."""
    with open(filepath, 'rb') as f:
        data = f.read()

    if data[0:10] != b'NeuroGuide':
        raise ValueError(f"Not a NeuroGuide file: {filepath}")

    strings = _scan_length_prefixed_strings(data, 10, 0x400)

    version = date = time_str = reference = ''
    time_end = None

    for offset, slen, text in strings:
        if '/' in text and len(text) == 10 and text[2] == '/':
            date = text
        elif ':' in text and len(text) == 8 and text[2] == ':':
            time_str = text
            time_end = offset + 4 + slen
        elif text in ('LinkEars', 'LinkedEars', 'Avg', 'Average',
                      'AvgRef', 'CzRef', 'A1Ref', 'A2Ref'):
            reference = text
        elif '.' in text and len(text) <= 10 and text[0].isdigit():
            version = text
        # NOTE: the first free-text string in a NeuroGuide header may be a
        # subject identifier. We deliberately do NOT read it — the output
        # EDF is labelled with the de-identified case ID instead.

    if time_end is None:
        raise ValueError(f"Cannot locate time field in header: {filepath}")

    sample_rate = struct.unpack_from('<H', data, time_end)[0]
    total_samples = struct.unpack_from('<I', data, time_end + 2)[0]
    n_channel_map = struct.unpack_from('<I', data, time_end + 6)[0]

    if sample_rate == 0 or total_samples == 0:
        raise ValueError(
            f"Invalid sample_rate={sample_rate} or total_samples={total_samples}")

    # Locate the data section: it starts right after a fixed 220-byte run of
    # 0xFF padding. Skipping the whole run would start one byte late when the
    # first data byte is itself 0xFF.
    search_from = time_end + 10 + n_channel_map * 2
    data_offset = None
    i = search_from
    while i < min(len(data), search_from + 8000):
        if data[i] == 0xFF:
            j = i
            while j < len(data) and data[j] == 0xFF:
                j += 1
            if j - i >= 10:
                data_offset = i + 220 if j - i >= 220 else j
                break
        i += 1

    if data_offset is None:
        raise ValueError(f"Cannot find data section in: {filepath}")

    if not reference:
        for kw in [b'LinkEars', b'LinkedEars', b'Avg', b'Average']:
            idx = data[:data_offset].find(kw)
            if idx >= 0:
                end_idx = data.index(b'\x00', idx)
                reference = data[idx:end_idx].decode('ascii', errors='replace')
                break

    raw = data[data_offset:]
    trim = len(raw) - (len(raw) % 4)
    all_floats = np.frombuffer(raw[:trim], dtype='<f4')

    # Count contiguous non-empty EEG channels from the start of the block.
    max_ch = len(all_floats) // total_samples
    n_eeg = 0
    for c in range(max_ch):
        chunk = all_floats[c * total_samples:(c + 1) * total_samples]
        if len(chunk) == total_samples and np.any(chunk != 0):
            if c == n_eeg:
                n_eeg += 1
            else:
                break
        else:
            break

    if n_eeg == 0:
        raise ValueError(f"No EEG data channels found in: {filepath}")

    channels = [all_floats[c * total_samples:(c + 1) * total_samples]
                for c in range(n_eeg)]

    ref_suffix = ''
    if reference == 'LinkEars':
        ref_suffix = '-LE'
    elif reference:
        ref_suffix = f'-{reference}'

    labels = []
    for c in range(n_eeg):
        base = NG_CHANNEL_NAMES_19[c] if c < len(NG_CHANNEL_NAMES_19) else f'Ch{c + 1}'
        labels.append(base + ref_suffix)

    # A channel stored as NaN (e.g. an electrode removed in NeuroGuide) cannot
    # be written to EDF. Set it to 0 and warn, so it is not written as garbage.
    for c in range(n_eeg):
        bad = ~np.isfinite(channels[c])
        if bad.any():
            channels[c] = np.where(bad, 0.0, channels[c]).astype(channels[c].dtype)
            print(f"  WARNING: {labels[c]} has {int(bad.sum())} non-finite samples "
                  f"(set to 0; mark it as a bad channel): {filepath}")

    return {
        'sample_rate': sample_rate,
        'total_samples': total_samples,
        'n_eeg': n_eeg,
        'channels': channels,
        'date': date,
        'time': time_str,
        'version': version,
        'reference': reference,
        'labels': labels,
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
            info = parse_ng(ng_path)
            duration = info['total_samples'] / info['sample_rate']
            print(f"  {case_id}: {info['n_eeg']}ch | {info['sample_rate']}Hz | "
                  f"{info['total_samples']} samples ({duration:.0f}s) | "
                  f"ref={info['reference'] or 'n/a'}")
            write_edf(edf_path, info, patient_id=case_id)
            print(f"    -> {case_id}.edf")
            converted += 1
        except Exception as e:                                    # noqa: BLE001
            print(f"  ERROR {case_id} ({ng_path}): {e}")
            failed += 1

    print(f"\nDone: {converted} converted, {failed} failed.")
    return 0 if failed == 0 else 2


if __name__ == '__main__':
    sys.exit(main())
