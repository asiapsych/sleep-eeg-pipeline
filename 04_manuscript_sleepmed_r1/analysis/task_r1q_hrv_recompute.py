"""從原始 ECG 重算 HRV；task_r1m 以這裡通過品管的結果為準。

NeuroKit2 清理、偵測 R 波，signal_fixpeaks 以 Kubios 方法校正異常心跳。
品管門檻在看 EEG–HRV 關係前訂定：模板相關中位數 ≥ .90、相關 ≥ .80 的
心跳至少 95%、校正比例 ≤ 5%、心率 40–120 bpm、RR 300–2000 ms。

EDF 的時間戳只能核對是否同次來訪，不能判斷 EEG 與 HRV 的記錄先後。
軟體值保留作對照；逐人檔不公開。"""

from __future__ import annotations
import csv
import datetime as dt
import io
import os
import sqlite3
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from _dynamic_n import DB_PATH, local_path
from _r1_common import AN, fmt_p, load_all, icc_a1 as icc_a1

warnings.filterwarnings("ignore")
RAW = Path(local_path("hrv_raw", "PSYML_HRV_RAW"))
EDF = Path(local_path("eeg_edf", "PSYML_EEG_EDF"))
sys.path.insert(0, str(RAW.parent / "scripts"))
OUT = AN / "S_R1q_hrv_recompute"
OUT.mkdir(parents=True, exist_ok=True)
FS = 2048
QC = dict(
    tpl_med=0.90, tpl_frac80=0.95, art_pct=5.0, hr=(40, 120), rr=(300, 2000)
)  # 訊號品管門檻（見檔頭第 5 點）
POWERLINE = 60


def edf_start(sid):
    p = EDF / f"{sid}.edf"
    if not p.exists():
        return None
    h = open(p, "rb").read(256)
    try:
        dd, mm, yy = h[168:176].decode("ascii").split(".")
        return dt.datetime(
            2000 + int(yy), int(mm), int(dd), *map(int, h[176:184].decode("ascii").split("."))
        )
    except Exception:
        return None


def session_info(folder):
    try:
        out = subprocess.check_output(
            ["mdb-export", str(folder / "ClSn.dat"), "Sessions"], stderr=subprocess.DEVNULL
        ).decode("utf-8", "replace")
        r = list(csv.DictReader(io.StringIO(out)))[0]
        return (
            dt.datetime.strptime(r["DateTime"], "%m/%d/%y %H:%M:%S"),
            r.get("FileName", ""),
            float(r.get("Duration") or 0),
        )
    except Exception:
        return None, "", np.nan


def recompute(sess):
    import neurokit2 as nk
    from read_session_dat import read_session_dat, read_statistics_dat

    s = read_session_dat(str(sess / "Session.dat"))
    cand = [c for c in s["channels"] if c["sample_rate"] == FS and c["sensor_type"] == 5]
    out = dict(duration=s["info"]["duration_sec"], has_ecg=bool(cand))
    if cand:
        try:
            x = np.asarray(cand[0]["data"], float)
            clean = nk.ecg_clean(x, sampling_rate=FS, method="neurokit", powerline=POWERLINE)
            _, info = nk.ecg_peaks(
                clean, sampling_rate=FS, method="neurokit", correct_artifacts=False
            )
            pk = np.asarray(info["ECG_R_Peaks"])
            art, fixed = nk.signal_fixpeaks(pk, sampling_rate=FS, iterative=True, method="Kubios")
            fixed = np.asarray(fixed)
            n_art = sum(len(art[k]) for k in ["ectopic", "missed", "extra", "longshort"])
            rr = np.diff(fixed) / FS * 1000
            a, b = int(0.25 * FS), int(0.40 * FS)
            keep = pk[(pk > a) & (pk < len(clean) - b)]
            beats = np.stack([clean[q - a : q + b] for q in keep])
            tpl = np.median(beats, 0)
            cc = np.array([np.corrcoef(bt, tpl)[0, 1] for bt in beats])
            t = nk.hrv_time({"ECG_R_Peaks": fixed}, sampling_rate=FS)
            f = nk.hrv_frequency({"ECG_R_Peaks": fixed}, sampling_rate=FS)
            out.update(
                n_beats=len(pk),
                art_pct=100 * n_art / len(pk),
                tpl_med=float(np.median(cc)),
                tpl_frac80=float(np.mean(cc >= 0.80)),
                rr_min=float(rr.min()),
                rr_max=float(rr.max()),
                HR=60000 / rr.mean(),
                SDNN=float(t["HRV_SDNN"][0]),
                RMSSD=float(t["HRV_RMSSD"][0]),
                pNN50=float(t["HRV_pNN50"][0]),
                LF=float(f["HRV_LF"][0]),
                HF=float(f["HRV_HF"][0]),
                LFHF=float(f["HRV_LFHF"][0]),
            )
        except Exception as e:
            out["error"] = str(e)[:80]
    st = sess / "Statistics.dat"
    if st.exists():
        try:
            out["software_HR"] = read_statistics_dat(str(st))["HR"]
        except Exception:
            pass
    return out


def qc_reasons(r):
    """回傳未通過品管的原因（空字串＝通過）。"""
    if not r.get("has_ecg", False) or pd.isna(r.get("RMSSD", np.nan)):
        return "no usable ECG channel"
    bad = []
    if r["tpl_med"] < QC["tpl_med"] or r["tpl_frac80"] < QC["tpl_frac80"]:
        bad.append("inconsistent beat morphology")
    if r["art_pct"] > QC["art_pct"]:
        bad.append("corrected beats > 5%")
    if not (QC["hr"][0] <= r["HR"] <= QC["hr"][1]):
        bad.append("mean heart rate outside 40–120 bpm")
    if r["rr_min"] < QC["rr"][0] or r["rr_max"] > QC["rr"][1]:
        bad.append("RR interval outside 300–2000 ms")
    return "; ".join(bad)


def main():
    df = load_all()
    with sqlite3.connect(DB_PATH) as con:
        h = pd.read_sql(
            "SELECT NEW_NUMBER AS subject, HRV_EKG_HR, HRV_SDNN_MS, HRV_RMSSD_MS, HRV_PNN50, HRV_HF, HRV_LF, PSQI_T1_C6 FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
            con,
        )
    for c in h.columns[1:]:
        h[c] = pd.to_numeric(h[c], errors="coerce")
    rows = []
    for sid in df["subject"]:
        folder = RAW / sid
        sess = folder / "Sess0001"
        if not (sess / "Session.dat").exists():
            continue
        r = dict(subject=sid)
        r.update(recompute(sess))
        when, prot, dur = session_info(folder)
        es = edf_start(sid)
        r.update(
            hrv_start=when,
            protocol=prot,
            eeg_start=es,
            gap_min=(when - es).total_seconds() / 60 if (when and es) else np.nan,
        )
        rows.append(r)
    rc = pd.DataFrame(rows)
    rc["qc_fail"] = [qc_reasons(r) for r in rc.to_dict("records")]
    rc["ok"] = rc["qc_fail"] == ""
    d = df.merge(rc, on="subject", how="left").merge(h, on="subject", how="left")
    d["has_raw"] = d["duration"].notna()
    d["ok"] = d["ok"].fillna(False).astype(bool)
    cols = [
        "subject",
        "has_raw",
        "ok",
        "qc_fail",
        "n_beats",
        "art_pct",
        "tpl_med",
        "tpl_frac80",
        "rr_min",
        "rr_max",
        "HR",
        "SDNN",
        "RMSSD",
        "pNN50",
        "LF",
        "HF",
        "LFHF",
        "software_HR",
        "HRV_EKG_HR",
        "HRV_SDNN_MS",
        "HRV_RMSSD_MS",
        "hrv_start",
        "eeg_start",
        "gap_min",
        "protocol",
    ]
    d[cols].to_csv(OUT / "_subject_level_internal.csv", index=False)

    # 品管結果（各組人數）
    qrows = []
    for g in ["HC", "Lower", "Higher", "All"]:
        x = d if g == "All" else d[d.k2 == g]
        raw = x[x.has_raw]
        row = dict(
            group=g,
            n_sample=len(x),
            n_with_raw_file=len(raw),
            n_pass=int(raw.ok.sum()),
            pct_pass_of_raw=100 * raw.ok.mean(),
        )
        for reason in [
            "no usable ECG channel",
            "inconsistent beat morphology",
            "corrected beats > 5%",
            "mean heart rate outside 40–120 bpm",
            "RR interval outside 300–2000 ms",
        ]:
            row[reason] = int(raw.qc_fail.fillna("").str.contains(reason, regex=False).sum())
        qrows.append(row)
    qc = pd.DataFrame(qrows)
    qc.to_csv(OUT / "01_qc.csv", index=False)
    from scipy.stats import chi2_contingency

    ct = pd.crosstab(d[d.has_raw].k2, d[d.has_raw].ok)
    chi_p = chi2_contingency(ct)[1] if ct.shape == (3, 2) else np.nan
    okd = d[d.ok]

    # 與軟體數值的對照（通過品管者）
    ag = []
    for new_, old, lab in [
        ("HR", "HRV_EKG_HR", "Heart rate, bpm"),
        ("SDNN", "HRV_SDNN_MS", "SDNN, ms"),
        ("RMSSD", "HRV_RMSSD_MS", "RMSSD, ms"),
    ]:
        x = okd[[new_, old]].dropna()
        x = x[x[old] > 0]
        a, b = x[old].to_numpy(), x[new_].to_numpy()
        df_ = b - a
        ag.append(
            dict(
                measure=lab,
                n=len(x),
                spearman_rho=spearmanr(a, b)[0],
                pearson_r_log=pearsonr(np.log(a), np.log(b))[0],
                median_software=float(np.median(a)),
                median_computed=float(np.median(b)),
                median_abs_diff=float(np.median(np.abs(df_))),
                n_within_10pct=int((np.abs(df_) <= 0.1 * a).sum()),
                n_diff_gt_10pct=int((np.abs(df_) > 0.1 * a).sum()),
            )
        )
    ag = pd.DataFrame(ag)
    ag.to_csv(OUT / "02_agreement.csv", index=False)
    hx = okd[["HR", "HRV_EKG_HR"]].dropna()
    n_hr_disc = int(((hx.HR - hx.HRV_EKG_HR).abs() > 3).sum())
    n_hr_double = int(((hx.HRV_EKG_HR / hx.HR) > 1.7).sum())
    fail = d[d.has_raw & ~d.ok]
    sw_fail = fail["HRV_RMSSD_MS"].dropna()

    # 時間
    g = d["gap_min"].dropna()
    tm = pd.DataFrame(
        [
            dict(
                n=len(g),
                same_day=int((g.abs() < 24 * 60).sum()),
                within_30_min=int((g.abs() <= 30).sum()),
                within_60_min=int((g.abs() <= 60).sum()),
                median_gap_min=g.median(),
                q25=g.quantile(0.25),
                q75=g.quantile(0.75),
                min=g.min(),
                max=g.max(),
            )
        ]
    )
    tm.to_csv(OUT / "03_timing.csv", index=False)

    L = [
        "# S_R1q — 由原始心電訊號計算 HRV（含訊號品管）",
        "",
        f"- 分析樣本 {len(df)} 人；有原始記錄檔 {int(d.has_raw.sum())} 人；通過訊號品管 {int(d.ok.sum())} 人。",
        f"- 記錄腳本：{rc.protocol.value_counts().head(3).to_dict()}；時長中位數 {rc.duration.median():.0f} 秒；取樣率 {FS} Hz。",
        f"- 品管門檻：心跳波形一致度的中位數 ≥ {QC['tpl_med']:.2f} 且相關 ≥ .80 的心跳 ≥ {QC['tpl_frac80']:.0%}；被校正的心跳 ≤ {QC['art_pct']:.0f}%；平均心率 {QC['hr'][0]}–{QC['hr'][1]} bpm；RR 間期 {QC['rr'][0]}–{QC['rr'][1]} ms。",
        "",
        "## 品管結果",
        "",
        "| 組別 | 樣本 | 有原始檔 | 通過 | 通過率 | 無可用心電 | 波形不一致 | 校正 > 5% | 心率超出範圍 | RR 超出範圍 |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in qc.to_dict("records"):
        L.append(
            f"| {r['group']} | {r['n_sample']} | {r['n_with_raw_file']} | {r['n_pass']} | {r['pct_pass_of_raw']:.1f}% | {r['no usable ECG channel']} | {r['inconsistent beat morphology']} | {r['corrected beats > 5%']} | {r['mean heart rate outside 40–120 bpm']} | {r['RR interval outside 300–2000 ms']} |"
        )
    L += [
        "",
        f"一筆記錄可能同時有多個原因。三組通過率的卡方檢定 p = {fmt_p(chi_p)}。",
        f"通過品管者：被校正的心跳中位數 {okd.art_pct.median():.2f}%（最大 {okd.art_pct.max():.1f}%）；心率 {okd.HR.min():.0f}–{okd.HR.max():.0f} bpm；RMSSD 中位數 {okd.RMSSD.median():.1f} ms（範圍 {okd.RMSSD.min():.1f}–{okd.RMSSD.max():.1f}）。",
        f"未通過品管的 {len(fail)} 筆中，{len(sw_fail)} 筆在總表仍有軟體算出的 RMSSD（中位數 {sw_fail.median():.1f} ms，最大 {sw_fail.max():.1f} ms）。",
        "",
        "## 與軟體數值的對照（通過品管者）",
        "",
        "| 指標 | n | Spearman rho | Pearson r（取對數後） | 軟體中位數 | 本腳本中位數 | 絕對差中位數 | 差在 10% 內的人數 | 差超過 10% 的人數 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in ag.itertuples():
        L.append(
            f"| {r.measure} | {r.n} | {r.spearman_rho:.2f} | {r.pearson_r_log:.2f} | {r.median_software:.1f} | {r.median_computed:.1f} | {r.median_abs_diff:.2f} | {r.n_within_10pct} | {r.n_diff_gt_10pct} |"
        )
    t0 = tm.iloc[0]
    L += [
        "",
        f"心率與軟體相差超過 3 bpm 的有 {n_hr_disc} 筆，其中軟體的心率約為本腳本兩倍的有 {n_hr_double} 筆。",
        "",
        "## HRV 與腦波記錄的起始時間差（HRV 減腦波，分鐘）",
        "",
        f"有時間資料 {int(t0.n)} 人；同一天 {int(t0.same_day)} 人；相差 30 分鐘內 {int(t0.within_30_min)} 人、60 分鐘內 {int(t0.within_60_min)} 人；中位數 {t0.median_gap_min:.1f}（四分位 {t0.q25:.1f} 至 {t0.q75:.1f}）。",
        "腦波檔的時間戳是建立檔案的時間，不能用來判斷先後；兩者實際上同時記錄。",
        "",
        "週期性 beta 與 HRV 的關係見 task_r1m_hrv_convergence.py（S_R1m_hrv）。",
    ]
    (OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
