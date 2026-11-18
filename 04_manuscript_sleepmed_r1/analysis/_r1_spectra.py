"""逐人、逐電極扣除 aperiodic 背景，保留 1–45 Hz 頻譜。

使用 S10 已擬合的 offset／exponent，不重新擬合 SpecParam。"""

from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.io import loadmat

from _r1_common import AN

PSD_DIR = AN / "S10_fooof_psd_input"
CH_RENAME = {"T7": "T3", "T8": "T4", "P7": "T5", "P8": "T6"}


def load_spectra(subjects):
    """回傳 freqs（1–45 Hz）、電極名、adj[人 × 電極 × 頻率]（去背景）、raw[人 × 電極 × 頻率]（log10 功率）。"""
    chp = pd.read_csv(AN / "S10_fooof_full" / "01_specparam_features.csv").set_index(
        ["subject", "channel"]
    )
    adj, raw, chn, freq = [], [], None, None
    for sid in subjects:
        m = loadmat(str(PSD_DIR / f"psd_{sid}.mat"), squeeze_me=True)
        f = np.asarray(m["freqs"], float)
        sel = (f >= 1.0) & (f <= 45.0)
        names = [
            CH_RENAME.get(str(c).upper(), str(c).upper()) for c in np.atleast_1d(m["ch_names"])
        ]
        lp = np.log10(np.asarray(m["psd"], float)[:, sel])
        ff = f[sel]
        ap = np.array([[chp.loc[(sid, c), "offset"], chp.loc[(sid, c), "exponent"]] for c in names])
        fl = lp - (ap[:, [0]] - ap[:, [1]] * np.log10(ff)[None, :])
        order = np.argsort(names)
        if chn is None:
            chn, freq = [names[i] for i in order], ff
        adj.append(fl[order])
        raw.append(lp[order])
    return freq, chn, np.asarray(adj), np.asarray(raw)


def band_mean(spec, freq, lo, hi):
    """spec[人 × 電極 × 頻率] 在 [lo, hi) 的平均，再對電極平均 → 每人一個值。"""
    sel = (freq >= lo) & (freq < hi) if hi < 45 else (freq >= lo) & (freq <= hi)
    return spec[:, :, sel].mean(2).mean(1)
