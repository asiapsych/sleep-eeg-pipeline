"""
Tests for ng2edf.py, using synthetic .ng files only (see ngsynth.py).

Run from the repository root:

    python -m unittest discover -s 01_ng2edf/tests -v
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import ng2edf  # noqa: E402
import ngsynth  # noqa: E402
from ngsynth import build_ng  # noqa: E402

N = 750                      # 2.5 s at 300 Hz (a multiple of 6, like real files)
LABELS_LE = [name + '-LE' for name in ng2edf.NG_CHANNEL_NAMES_19]

# SHA-256 of the EDF that the previous converter (origin/main e697e68)
# writes for these synthetic files. The new reader must reproduce them.
GOLDEN_EDF_SHA256 = {
    'default': 'ca0796598276fbd227deda2b0aeff9f00a6c78df04eefa0693da93deb22a8543',
    'first_byte_ff': 'ef0570ef425e48d9005a7803189111f3668190ffeb585df86c1e466582ec6162',
    'nan_pz': '0b8a66e5ddc33a9620f1bf9173e76fe13aa800f228486654462201556121a8ac',
}


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='ng2edf_test_')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write_ng(self, buf, case='S001', name='rec.ng'):
        folder = os.path.join(self.tmp, 'in', case)
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, name)
        with open(path, 'wb') as f:
            f.write(buf)
        return path

    def parse(self, buf, **kw):
        path = self.write_ng(buf)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            info = ng2edf.parse_ng(path, **kw)
        return info, out.getvalue()

    def run_cli(self, *extra):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = ng2edf.main(['-i', os.path.join(self.tmp, 'in'),
                              '-o', os.path.join(self.tmp, 'out'), *extra])
        return rc, out.getvalue()

    def out_path(self, name):
        return os.path.join(self.tmp, 'out', name)


def read_edf(path):
    """Minimal EDF reader: (header fields, labels, physical data)."""
    with open(path, 'rb') as f:
        raw = f.read()
    ns = int(raw[252:256])
    hdr = {
        'version': raw[0:8].decode().strip(),
        'patient': raw[8:88].decode().strip(),
        'recording': raw[88:168].decode().strip(),
        'startdate': raw[168:176].decode(),
        'starttime': raw[176:184].decode(),
        'header_bytes': int(raw[184:192]),
        'n_records': int(raw[236:244]),
        'duration': raw[244:252].decode().strip(),
        'ns': ns,
    }

    def field(start, width):
        base = 256 + start * ns
        return [raw[base + i * width:base + (i + 1) * width].decode().strip()
                for i in range(ns)]

    labels = field(0, 16)
    pmin = np.array(field(16 + 80 + 8, 8), dtype=float)
    pmax = np.array(field(16 + 80 + 16, 8), dtype=float)
    dmin = np.array(field(16 + 80 + 24, 8), dtype=float)
    dmax = np.array(field(16 + 80 + 32, 8), dtype=float)
    spr = int(field(16 + 80 + 40 + 80, 8)[0])
    dig = np.frombuffer(raw[hdr['header_bytes']:], dtype='<i2')
    dig = dig.reshape(hdr['n_records'], ns, spr).transpose(1, 0, 2).reshape(ns, -1)
    gain = (pmax - pmin) / (dmax - dmin)
    phys = (dig - dmin[:, None]) * gain[:, None] + pmin[:, None]
    return hdr, labels, phys, gain[0]


# ── header / section A ───────────────────────────────────────────────────

class HeaderTests(TempDirTestCase):
    def test_standard_layout_offsets(self):
        buf = build_ng(n=N)
        hdr = ng2edf.parse_header(buf)
        # age 5 chars, date 10, label 3, time 8, all else empty -> 569 bytes
        self.assertEqual(hdr['data_offset'], 569)
        self.assertEqual(hdr['tail_offset'], 569 + 4 * 21 * N + 21)
        self.assertEqual(len(buf), hdr['tail_offset'] + 5 * N)
        self.assertEqual(hdr['subject_offset'], 33)
        self.assertEqual((hdr['sample_rate'], hdr['total_samples'], hdr['n_channels']),
                         (300, N, 21))
        self.assertEqual(hdr['channel_codes'], list(range(19)) + [55, 56])
        self.assertEqual(hdr['ref_code'], 55)
        self.assertEqual(hdr['reference'], 'LinkEars')
        self.assertEqual((hdr['date'], hdr['time'], hdr['age']),
                         ('01/31/2000', '12:34:56', '34.00'))
        self.assertFalse(hdr['section_a_deviates'])
        self.assertFalse(hdr['section_a_ambiguous'])

    def test_nonempty_string_without_flag_byte(self):
        text = b'SYNTHETIC FREE TEXT'
        strings = {**ngsynth.DEFAULT_STRINGS, 1: text}
        buf = build_ng(n=N, strings=strings)
        info, out = self.parse(buf)
        self.assertEqual(info['data_offset'], 569 + len(text))
        self.assertEqual((info['date'], info['time'], info['age']),
                         ('01/31/2000', '12:34:56', '34.00'))
        self.assertEqual(info['labels'], LABELS_LE)
        self.assertEqual(out, '')
        np.testing.assert_array_equal(np.vstack(info['channels']), ngsynth.synthetic_eeg(N))

    def test_legacy_flag_before_nonempty_string(self):
        # Older reading of section A: a 0x01 byte before a non-empty string.
        strings = {**ngsynth.DEFAULT_STRINGS, 4: b'XYZ'}
        buf = build_ng(n=N, strings=strings, legacy_flags={4})
        info, out = self.parse(buf)
        self.assertIn('located by search', out)
        self.assertEqual((info['date'], info['time']), ('01/31/2000', '12:34:56'))
        np.testing.assert_array_equal(np.vstack(info['channels']), ngsynth.synthetic_eeg(N))

    def test_empty_strings_with_length_zero(self):
        buf = build_ng(n=N, len0_slots={1, 2, 8})
        hdr = ng2edf.parse_header(buf)
        self.assertFalse(hdr['section_a_deviates'])
        self.assertEqual(hdr['data_offset'], 569 - 3)
        self.assertEqual(hdr['date'], '01/31/2000')

    def test_field_byte_before_label_any_value(self):
        for value in (0, 2, 0xFF):
            hdr = ng2edf.parse_header(build_ng(n=N, field_byte=value))
            self.assertFalse(hdr['section_a_deviates'], value)
            self.assertEqual(hdr['time'], '12:34:56')

    def test_field_byte_missing(self):
        hdr = ng2edf.parse_header(build_ng(n=N, field_byte=None))
        self.assertTrue(hdr['section_a_deviates'])
        self.assertEqual((hdr['date'], hdr['time']), ('01/31/2000', '12:34:56'))

    def test_longer_header_uses_pointer(self):
        strings = {**ngsynth.DEFAULT_STRINGS, 12: b'LABEL-0000000123'}
        info, _ = self.parse(build_ng(n=N, strings=strings))
        self.assertEqual(info['data_offset'], 569 + 13)
        np.testing.assert_array_equal(np.vstack(info['channels']), ngsynth.synthetic_eeg(N))

    def test_malformed_date_and_time_are_not_used(self):
        strings = {**ngsynth.DEFAULT_STRINGS, 7: b'2000-01-31', 13: b'noon'}
        info, out = self.parse(build_ng(n=N, strings=strings))
        self.assertEqual((info['date'], info['time']), ('', ''))
        self.assertIn('date of test is not MM/DD/YYYY', out)
        self.assertIn('time of test is not HH:MM:SS', out)
        self.assertNotIn('2000-01-31', out)

    def test_bad_magic(self):
        buf = b'NotNeuroGd' + build_ng(n=N)[10:]
        with self.assertRaises(ng2edf.NGFormatError):
            self.parse(buf)

    def test_section_pointer_mismatch(self):
        buf = build_ng(n=N)
        off_a, off_b, _, _ = ngsynth.offsets(buf)
        bad = build_ng(n=N, pointer_patch={0x15: off_b + 1})
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'section A'):
            ng2edf.parse_header(bad)
        bad = build_ng(n=N, pointer_patch={0x19: ngsynth.offsets(buf)[2] - 2})
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'section B ends'):
            ng2edf.parse_header(bad)

    def test_unrecognised_preamble_constants_warn(self):
        _, out = self.parse(build_ng(n=N, consts=(3, 4, 8)))
        self.assertIn('unrecognised preamble constants', out)


# ── data section and channel table ───────────────────────────────────────

class DataTests(TempDirTestCase):
    def test_channels_labels_and_average_check(self):
        info, out = self.parse(build_ng(n=N))
        self.assertEqual(out, '')
        self.assertEqual(info['n_eeg'], 19)
        self.assertEqual(info['labels'], LABELS_LE)
        self.assertEqual(info['channel_codes'], list(range(19)))
        self.assertEqual(info['checks']['average'], 'ok')
        self.assertNotIn('version', info)
        np.testing.assert_array_equal(np.vstack(info['channels']), ngsynth.synthetic_eeg(N))

    def test_first_data_byte_ff(self):
        eeg = ngsynth.synthetic_eeg(N)
        eeg[0, 0] = ngsynth.FIRST_BYTE_FF
        buf = build_ng(n=N, eeg=eeg)
        off_data = ngsynth.offsets(buf)[2]
        self.assertEqual(buf[off_data - 220:off_data + 1], b'\xff' * 221)
        info, _ = self.parse(buf)
        self.assertEqual(info['channels'][0][0], ngsynth.FIRST_BYTE_FF)
        np.testing.assert_array_equal(np.vstack(info['channels']), eeg)

    def test_all_nan_channel_set_to_zero_and_warned(self):
        eeg = ngsynth.synthetic_eeg(N)
        eeg[18, :] = ngsynth.NAN32
        info, out = self.parse(build_ng(n=N, eeg=eeg))
        self.assertIn(f'WARNING: Pz-LE has {N} non-finite samples '
                      f'(set to 0; mark it as a bad channel)', out)
        self.assertTrue(np.all(info['channels'][18] == 0))
        self.assertTrue(info['checks']['average'].startswith('skipped'))
        np.testing.assert_array_equal(np.vstack(info['channels'][:18]), eeg[:18])

    def test_average_mismatch_is_an_error(self):
        avg = ngsynth.average_block(ngsynth.synthetic_eeg(N))
        avg[100] += 1.0
        buf = build_ng(n=N, blocks={56: avg})
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'average channel'):
            self.parse(buf)
        info, out = self.parse(buf, lenient=True)
        self.assertEqual(info['checks']['average'], 'mismatch')
        self.assertIn('average channel', out)

    def test_misread_channel_table_is_caught(self):
        # The data blocks are shifted by one against the channel table (the
        # zero block first): the average check must refuse the file.
        eeg = ngsynth.synthetic_eeg(N)
        blocks = {0: np.zeros(N)}
        blocks.update({code: eeg[code - 1] for code in range(1, 19)})
        blocks[55] = eeg[18]
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'average channel'):
            self.parse(build_ng(n=N, blocks=blocks))

    def test_channels_follow_the_codes(self):
        codes = [1, 0] + list(range(2, 19)) + [55, 56]
        info, _ = self.parse(build_ng(n=N, codes=codes))
        self.assertEqual(info['labels'][:3], ['FP2-LE', 'FP1-LE', 'F3-LE'])
        eeg = ngsynth.synthetic_eeg(N)
        np.testing.assert_array_equal(info['channels'][0], eeg[1])
        self.assertEqual(info['checks']['average'], 'ok')

    def test_unknown_codes_are_ignored(self):
        codes = list(range(19)) + [55, 60, 56]
        info, out = self.parse(build_ng(n=N, codes=codes))
        self.assertEqual(info['n_eeg'], 19)
        self.assertIn('unknown electrode codes [60]', out)

    def test_fewer_eeg_channels(self):
        codes = list(range(16)) + [55, 56]
        eeg = ngsynth.synthetic_eeg(N)
        avg = ngsynth.average_block(eeg[:16])
        info, _ = self.parse(build_ng(n=N, codes=codes, blocks={56: avg}))
        self.assertEqual(info['n_eeg'], 16)
        self.assertEqual(info['labels'][-1], 'T6-LE')
        self.assertEqual(info['checks']['average'], 'ok')

    def test_no_eeg_codes(self):
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'No EEG data channels'):
            self.parse(build_ng(n=N, codes=[55, 56]))

    def test_unknown_montage_name(self):
        info, out = self.parse(build_ng(n=N, montage=b'MyMontage'))
        self.assertEqual(info['labels'][0], 'FP1')
        self.assertEqual(info['reference'], '')
        self.assertIn('montage name not recognised', out)
        self.assertNotIn('MyMontage', out)

    def test_average_montage_suffix(self):
        info, _ = self.parse(build_ng(n=N, montage=b'Avg'))
        self.assertEqual(info['labels'][0], 'FP1-Avg')


# ── file size and section offsets ────────────────────────────────────────

class SizeTests(TempDirTestCase):
    def test_truncated_tail(self):
        buf = build_ng(n=N, events=[(10, 2)], truncate=7)
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'file size'):
            self.parse(buf)
        info, out = self.parse(buf, lenient=True)
        self.assertIn('file size', out)
        self.assertEqual([e['record'] for e in info['events']], [10])

    def test_extra_trailing_bytes(self):
        buf = build_ng(n=N, extra=b'\x00' * 15)
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'expected'):
            self.parse(buf)
        info, _ = self.parse(buf, lenient=True)
        np.testing.assert_array_equal(np.vstack(info['channels']), ngsynth.synthetic_eeg(N))

    def test_truncated_data_is_always_an_error(self):
        buf = build_ng(n=N, truncate=5 * N + 21 + 100)
        for lenient in (False, True):
            with self.assertRaises(ng2edf.NGFormatError):
                self.parse(buf, lenient=lenient)

    def test_tail_pointer_inside_data(self):
        buf = build_ng(n=N)
        _, _, off_data, off_tail = ngsynth.offsets(buf)
        bad = build_ng(n=N, pointer_patch={0x1d: off_tail - 100})
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'past the tail offset'):
            self.parse(bad, lenient=True)

    def test_pointers_out_of_order(self):
        bad = build_ng(n=N, pointer_patch={0x19: 10})
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'out of order'):
            ng2edf.parse_header(bad)

    def test_nonzero_gap_only_warns(self):
        info, out = self.parse(build_ng(n=N, gap=b'\x00' * 20 + b'\x01'))
        self.assertIn('bytes between data and tail', out)
        self.assertEqual(info['n_eeg'], 19)

    def test_zero_samples(self):
        with self.assertRaisesRegex(ng2edf.NGFormatError, 'total_samples=0'):
            ng2edf.parse_header(build_ng(n=0, eeg=np.zeros((19, 0))))


# ── tail section / events ────────────────────────────────────────────────

class EventTests(TempDirTestCase):
    def test_events_at_known_tail_offsets(self):
        last = 5 * N // 6 - 1
        buf = build_ng(n=N, events=[(0, 2), (123, 2), (last, 7)])
        info, out = self.parse(buf)
        self.assertEqual(out, '')
        tail_offset = ngsynth.offsets(buf)[3]
        got = [(e['record'], e['tail_offset'], e['file_offset'], e['code'], e['kind'])
               for e in info['events']]
        self.assertEqual(got, [
            (0, 4, tail_offset + 4, 2, 'event'),
            (123, 742, tail_offset + 742, 2, 'event'),
            (last, 6 * last + 4, tail_offset + 6 * last + 4, 7, 'event'),
        ])
        for e in info['events']:
            self.assertEqual(buf[e['file_offset']], e['code'])
            self.assertNotIn('time', ' '.join(e))

    def test_event_time_needs_an_explicit_base(self):
        ev = ng2edf.parse_events(build_ng(n=N, events=[(500, 2)])[-5 * N:])[0]
        # 300 Hz: 250, 200 and 300 records per second
        self.assertEqual(ng2edf.event_time_s(ev, 300, 'span'), 2.0)
        self.assertEqual(ng2edf.event_time_s(ev, 300, 'legacy'), 2.5)
        self.assertAlmostEqual(ng2edf.event_time_s(ev, 300, 'sample'), 500 / 300)
        # the candidates scale with the sample rate
        self.assertEqual(ng2edf.event_time_s(ev, 250, 'sample'), 2.0)
        with self.assertRaises(TypeError):
            ng2edf.event_time_s(ev, 300)     # pylint: disable=no-value-for-parameter
        with self.assertRaises(ValueError):
            ng2edf.event_time_s(ev, 300, 250)

    def test_unrecognised_tail_record(self):
        info, out = self.parse(build_ng(n=N, tail_patch={61: 9}))
        self.assertEqual(info['events'][0]['kind'], 'unknown')
        self.assertEqual(info['events'][0]['record'], 10)
        self.assertEqual(info['events'][0]['raw_hex'], '000900000000')
        self.assertIn('1 tail record(s) with unrecognised content', out)

    def test_no_events(self):
        info, _ = self.parse(build_ng(n=N))
        self.assertEqual(info['events'], [])

    def test_events_tsv(self):
        self.write_ng(build_ng(n=N, events=[(300, 2), (601, 2)]))
        rc, out = self.run_cli('--events')
        self.assertEqual(rc, 0)
        self.assertIn('-> S001_events.tsv (2 record(s))', out)
        with open(self.out_path('S001_events.tsv'), newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(list(rows[0]), ng2edf.EVENT_TSV_COLUMNS)
        self.assertEqual([r['record'] for r in rows], ['300', '601'])
        self.assertEqual(rows[0]['tail_offset'], str(6 * 300 + 4))
        self.assertEqual(rows[0]['time_s_span'], '1.200000')
        self.assertEqual(rows[0]['time_s_legacy'], '1.500000')
        self.assertEqual(rows[0]['time_s_sample'], '1.000000')
        self.assertEqual(rows[1]['time_s_sample'], '2.003333')
        self.assertEqual(rows[0]['kind'], 'event')

    def test_events_tsv_header_only_when_empty(self):
        self.write_ng(build_ng(n=N))
        self.run_cli('--events')
        with open(self.out_path('S001_events.tsv')) as f:
            self.assertEqual(f.read(), '\t'.join(ng2edf.EVENT_TSV_COLUMNS) + '\n')

    def test_no_tsv_without_flag(self):
        self.write_ng(build_ng(n=N, events=[(300, 2)]))
        rc, out = self.run_cli()
        self.assertEqual(rc, 0)
        self.assertEqual(os.listdir(os.path.join(self.tmp, 'out')), ['S001.edf'])
        self.assertNotIn('events', out)


# ── EDF output and CLI ───────────────────────────────────────────────────

class EdfTests(TempDirTestCase):
    def convert(self, buf, *extra):
        self.write_ng(buf)
        rc, out = self.run_cli(*extra)
        self.assertEqual(rc, 0, out)
        return self.out_path('S001.edf'), out

    def test_golden_bytes_match_previous_converter(self):
        eeg_ff = ngsynth.synthetic_eeg(N)
        eeg_ff[0, 0] = ngsynth.FIRST_BYTE_FF
        eeg_nan = ngsynth.synthetic_eeg(N)
        eeg_nan[18, :] = ngsynth.NAN32
        cases = [('default', 'S001', build_ng(n=N)),
                 ('first_byte_ff', 'S002', build_ng(n=N, eeg=eeg_ff)),
                 ('nan_pz', 'S003', build_ng(n=N, eeg=eeg_nan))]
        for _, case, buf in cases:
            self.write_ng(buf, case=case)
        rc, out = self.run_cli()
        self.assertEqual(rc, 0, out)
        for name, case, _ in cases:
            with self.subTest(name):
                with open(self.out_path(f'{case}.edf'), 'rb') as f:
                    self.assertEqual(hashlib.sha256(f.read()).hexdigest(),
                                     GOLDEN_EDF_SHA256[name])

    def test_edf_header_and_data(self):
        path, out = self.convert(build_ng(n=N))
        self.assertIn('S001: 19ch | 300Hz | 750 samples (2s) | ref=LinkEars', out)
        hdr, labels, phys, gain = read_edf(path)
        self.assertEqual(hdr['version'], '0')
        self.assertEqual(hdr['patient'], 'S001')
        self.assertEqual(hdr['recording'], '')
        self.assertEqual(hdr['startdate'], '31.01.00')
        self.assertEqual(hdr['starttime'], '12.34.56')
        self.assertEqual(hdr['n_records'], 3)
        self.assertEqual(hdr['ns'], 19)
        self.assertEqual(labels, LABELS_LE)
        err = np.abs(phys[:, :N] - ngsynth.synthetic_eeg(N))
        self.assertLessEqual(err.max(), gain / 2 + 1e-9)

    def test_header_text_never_reaches_outputs(self):
        secret = b'SYNTHETIC-PERSON-NAME'
        strings = {**ngsynth.DEFAULT_STRINGS, 1: secret, 11: b'synthetic comment'}
        path, out = self.convert(build_ng(n=N, strings=strings, events=[(5, 2)]),
                                 '--events')
        with open(path, 'rb') as f:
            edf = f.read()
        with open(self.out_path('S001_events.tsv'), 'rb') as f:
            tsv = f.read()
        for needle in (secret, b'synthetic comment', b'999', b'34.00'):
            self.assertNotIn(needle, edf[:256 + 19 * 256])
            self.assertNotIn(needle, tsv)
            self.assertNotIn(needle.decode(), out)

    def test_dry_run_writes_nothing(self):
        self.write_ng(build_ng(n=N))
        rc, out = self.run_cli('--dry-run', '--events')
        self.assertEqual(rc, 0)
        self.assertIn('S001.edf  <-  ', out)
        self.assertEqual(os.listdir(os.path.join(self.tmp, 'out')), [])

    def test_failed_file_does_not_stop_the_batch(self):
        self.write_ng(build_ng(n=N), case='S001')
        self.write_ng(build_ng(n=N, truncate=3), case='S002')
        rc, out = self.run_cli()
        self.assertEqual(rc, 2)
        self.assertIn('ERROR S002', out)
        self.assertIn('Done: 1 converted, 1 failed.', out)
        rc, out = self.run_cli('--lenient')
        self.assertEqual(rc, 0)
        self.assertIn('Done: 2 converted, 0 failed.', out)


if __name__ == '__main__':
    unittest.main()
