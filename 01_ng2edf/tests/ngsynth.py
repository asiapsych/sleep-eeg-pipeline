"""
Builder for SYNTHETIC NeuroGuide ``.ng`` files, used by the tests.

The files follow the layout in ``../NG_FORMAT.md``. Every value (signals,
subject strings, events) is made up; no real recording is involved. Files
are only ever written to temporary directories.
"""

from __future__ import annotations

import struct

import numpy as np

EEG_CODES = list(range(19))
CODE_LE, CODE_AVG = 55, 56
DEFAULT_CODES = EEG_CODES + [CODE_LE, CODE_AVG]

# Made-up subject information (slot number -> bytes). Slot 3 = age,
# 7 = date of test, 12 = label, 13 = time of test; the rest stay empty.
DEFAULT_STRINGS = {3: b'34.00', 7: b'01/31/2000', 12: b'999', 13: b'12:34:56'}

NAN_BITS = 0xFFC00000                     # NaN pattern NeuroGuide writes
NAN32 = np.array([NAN_BITS], dtype='<u4').view('<f4')[0]
FIRST_BYTE_FF = np.frombuffer(b'\xff\x00\x80\x3f', dtype='<f4')[0]  # ~1.00003


def ng_string(text, empty_len0=False):
    """u32 length (incl. NUL) + bytes; empty = length 1 + NUL (or length 0)."""
    if not text:
        return struct.pack('<I', 0) if empty_len0 else struct.pack('<IB', 1, 0)
    return struct.pack('<I', len(text) + 1) + text + b'\x00'


def synthetic_eeg(n, n_eeg=19):
    """Deterministic test signal; every value is a multiple of 1/8 uV, so it
    is exact in float32 and the EDF output does not depend on the platform."""
    i = np.arange(n, dtype=np.int64)
    rows = [((i * (7919 + 2 * c) + 104729 * (c + 1)) % 2001 - 1000) / 8.0
            for c in range(n_eeg)]
    return np.array(rows, dtype=np.float32)


def average_block(eeg):
    """float32(sum in double / n_eeg): NeuroGuide's code-56 block."""
    return (eeg.astype(np.float64).sum(axis=0) / eeg.shape[0]).astype(np.float32)


def build_ng(n=750, fs=300, eeg=None, codes=None, blocks=None, strings=None,
             field_byte=1, legacy_flags=(), len0_slots=(), consts=(2, 4, 8),
             montage=b'LinkEars', ref_code=CODE_LE, t3_count=110,
             gap=b'\x00' * 21, events=(), tail_patch=None, tail=None,
             extra=b'', truncate=0, pointer_patch=None):
    """Return the bytes of a synthetic ``.ng`` file.

    eeg           19 x n float32 array indexed by electrode code (default:
                  :func:`synthetic_eeg`)
    codes         electrode code of every data block (default 0..18, 55, 56)
    blocks        {code: array} to override individual data blocks
    strings       {slot: bytes} for section A (default DEFAULT_STRINGS)
    field_byte    value of the u8 before string #12, or None to omit it
    legacy_flags  slots preceded by a 0x01 "flag" byte (older reading)
    len0_slots    empty slots written as length 0 instead of length 1 + NUL
    events        [(record, code)] -> tail byte 6*record + 4 = code
    tail_patch    {tail offset: byte value}
    tail          explicit tail bytes (default 5*n zero bytes + events)
    extra         bytes appended after the tail
    truncate      number of bytes cut from the end
    pointer_patch {0x11|0x15|0x19|0x1d: value} to corrupt a section offset
    """
    codes = list(DEFAULT_CODES if codes is None else codes)
    eeg = synthetic_eeg(n) if eeg is None else np.asarray(eeg, dtype=np.float32)
    strings = dict(DEFAULT_STRINGS if strings is None else strings)

    sec_a = b''
    for slot in range(1, 14):
        if slot == 12 and field_byte is not None:
            sec_a += bytes([field_byte])
        if slot in legacy_flags:
            sec_a += b'\x01'
        sec_a += ng_string(strings.get(slot, b''),
                           empty_len0=slot in legacy_flags or slot in len0_slots)

    sec_b = struct.pack('<HIH', fs, n, len(codes))
    sec_b += struct.pack(f'<{len(codes)}H', *codes)
    sec_b += struct.pack('<HB', ref_code, 0x47)
    for _ in range(2):
        sec_b += struct.pack('<H', 19)
        sec_b += b''.join(struct.pack('<HH', i, ref_code) for i in range(19))
    sec_b += ng_string(montage)
    sec_b += struct.pack('<H', t3_count) + struct.pack('<h', -1) * t3_count

    data = b''
    for code in codes:
        if blocks is not None and code in blocks:
            blk = np.asarray(blocks[code], dtype='<f4')
        elif code < 19:
            blk = eeg[code]
        elif code == CODE_AVG:
            blk = average_block(eeg)
        else:
            blk = np.zeros(n, dtype='<f4')
        data += blk.astype('<f4').tobytes()

    off_a = 33
    off_b = off_a + len(sec_a)
    off_data = off_b + len(sec_b)
    off_tail = off_data + len(data) + len(gap)
    pointers = {0x11: off_a, 0x15: off_b, 0x19: off_data, 0x1d: off_tail}
    pointers.update(pointer_patch or {})
    pre = b'NeuroGuide\x00' + struct.pack('<3H', *consts)
    pre += struct.pack('<4I', pointers[0x11], pointers[0x15],
                       pointers[0x19], pointers[0x1d])
    assert len(pre) == 33

    if tail is None:
        t = bytearray(5 * n)
        for record, code in events:
            t[6 * record + 4] = code
        for off, val in (tail_patch or {}).items():
            t[off] = val
        tail = bytes(t)

    out = pre + sec_a + sec_b + data + gap + tail + extra
    return out[:len(out) - truncate] if truncate else out


def offsets(buf):
    """(section A, section B, data, tail) offsets stored in a file."""
    return struct.unpack_from('<4I', buf, 0x11)
