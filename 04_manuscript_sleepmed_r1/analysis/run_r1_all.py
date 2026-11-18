"""依序重跑 R1；每一步的 log 與雜湊記在 work/runs/<run_id>。

直接執行會先同步前處理頻譜，再跑分析及圖表。
  --no-sync                    沿用目前的頻譜檔
  --from task_r1a_severity_form.py   從指定腳本接著跑

一般步驟失敗即停止；LOCAL 步驟失敗會記錄後略過。"""

from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path

BASE = Path(__file__).resolve().parent
AN = BASE / "outputs" / "補充資料" / "analysis"
from _dynamic_n import local_path

PIPE = Path(
    local_path("eeg_pipeline_root", "EEG_PIPELINE_ROOT")
)
PSD_IN = AN / "S10_fooof_psd_input"
UPSTREAM = [
    "task_b_fooof_full.py",
    "task_f_recompute_relpower.py",
    "task_e_bootstrap_stability.py",
    "task_e2_bootstrap_selectK.py",
]
MAIN = [
    "task_r1a_severity_form.py",
    "task_r1b_k234_comparison.py",
    "task_r1c_sensitivity_emg_vigilance.py",
    "task_r1d_preproc_qc.py",
    "task_r1e_clustering_robustness.py",
    "task_r1f_pc_biomarkers.py",
    "task_r1g_followup.py",
    "task_r1h_beta_characterisation.py",
    "task_r1i_mood_adjustment.py",
    "task_r1k_recruitment_flags.py",
    "task_r1q_hrv_recompute.py",
    "task_r1m_hrv_convergence.py",
    "task_r1n_beta_subbands.py",
    "task_r1o_sleep_medication.py",
    "task_r1s_vigilance_epochs.py",
    "task_r1u_shape_robust.py",
    "task_r1v_sleep_onset_check.py",
    "task_r1w_pc_interpretation.py",
]
LOCAL = [
    "task_r1j_manual_vs_auto.py",
    "task_r1p_picard_comparison.py",
    "task_r1t_surface_laplacian.py",
]  # 需要本機的原始資料
SUBMISSION = [
    "sub_table1.py",
    "sub_table2.py",
    "sub_fig1.py",
    "sub_fig2.py",
    "sub_fig3.py",
    "sub_fig4.py",
    "sub_fig5.py",
    "sub_fig6.py",
    "sub_figS1.py",
    "sub_figS2.py",
    "sub_tableS1.py",
    "sub_tableS2.py",
    "sub_tableS3.py",
    "sub_tableS4.py",
    "sub_tableS5.py",
    "sub_response_tables.py",
]  # 送件用圖表，一張一支
FINAL = (
    ["make_r1_report_figures.py", "make_r1_report_tables.py"]
    + SUBMISSION
    + ["make_r1_ssot.py", "sub_sync_notes.py"]
)
PKGS = [
    "numpy",
    "pandas",
    "scipy",
    "statsmodels",
    "scikit-learn",
    "stepmix",
    "fooof",
    "neurokit2",
    "mne",
    "matplotlib",
    "openpyxl",
]


def sha(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def sha_dir(d, pattern):
    h = hashlib.sha256()
    files = sorted(Path(d).glob(pattern))
    for p in files:
        h.update(p.name.encode())
        h.update(sha(p, 64).encode())
    return h.hexdigest()[:16], len(files)


def sync_psd(man):
    src = PIPE / "outputs" / "psd_mat"
    files = sorted(src.glob("psd_*.mat"))
    if not files:
        raise SystemExit(f"找不到前處理輸出：{src}")
    if PSD_IN.exists():
        shutil.rmtree(PSD_IN)
    PSD_IN.mkdir(parents=True)
    for p in files:
        shutil.copy2(p, PSD_IN / p.name)
    man["psd_sync"] = dict(source=str(src), n_files=len(files))


def pipeline_stamp():
    mdir = PIPE / "outputs" / "manifests"
    out = dict(root=str(PIPE))
    try:
        p = sorted(mdir.glob("*.json"))[0]
        m = json.loads(p.read_text())
        vs = m.get("version_stamp", {})
        out.update(
            {
                k: vs.get(k)
                for k in (
                    "pipeline_version",
                    "git_hash",
                    "matlab_version",
                    "eeglab_version",
                    "timestamp",
                )
            }
        )
        out["plugins"] = vs.get("plugins")
        s2 = m.get("step2", {})
        out["ica_algorithm"] = s2.get("ica_algorithm")
        out["n_manifests"] = len(list(mdir.glob("*.json")))
    except Exception as e:
        out["error"] = str(e)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", default=None)
    ap.add_argument("--no-sync", action="store_true")
    args = ap.parse_args()
    run_id = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    rd = BASE / "work" / "runs" / run_id
    rd.mkdir(parents=True)
    man = dict(
        run_id=run_id,
        started=dt.datetime.now().isoformat(timespec="seconds"),
        host=platform.node(),
        platform=platform.platform(),
        python=sys.version.split()[0],
        sample_env=os.environ.get("PSYML_SAMPLE", "audited"),
        steps=[],
    )
    man["packages"] = {}
    for k in PKGS:
        try:
            man["packages"][k] = metadata.version(k)
        except Exception:
            man["packages"][k] = None
    try:
        man["git_commit"] = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=BASE, stderr=subprocess.DEVNULL
            )
            .decode()
            .strip()
        )
    except Exception:
        man["git_commit"] = None
    steps = UPSTREAM + MAIN + LOCAL + FINAL
    if args.start:
        steps = steps[steps.index(args.start) :]
    if not args.no_sync and not args.start:
        sync_psd(man)
    man["pipeline"] = pipeline_stamp()
    h, n = sha_dir(PSD_IN, "psd_*.mat")
    man["inputs"] = dict(psd_dir=str(PSD_IN), psd_n=n, psd_sha=h)
    sys.path.insert(0, str(BASE))
    from _dynamic_n import DB_PATH, excluded_ids

    man["inputs"].update(
        db=str(DB_PATH),
        db_sha=sha(DB_PATH),
        db_mtime=dt.datetime.fromtimestamp(os.path.getmtime(DB_PATH)).isoformat(timespec="seconds"),
    )
    man["sample_definition"] = dict(n_excluded=len(excluded_ids()), excluded=excluded_ids())
    man["scripts"] = {p.name: sha(p) for p in sorted(BASE.glob("*.py"))}

    def save():
        (rd / "run_manifest.json").write_text(
            json.dumps(man, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    save()
    ok = True
    for s in steps:
        if not (BASE / s).exists():          # 這份程式沒有附的步驟（例如只在本機用的腳本）直接略過
            man["steps"].append(dict(script=s, status="skipped (not in this copy)", seconds=0, log=""))
            print(f"{s}: skipped (not in this copy)", flush=True)
            save()
            continue
        t0 = time.time()
        logp = rd / (s.replace(".py", "") + ".log")
        with open(logp, "w", encoding="utf-8") as fh:
            r = subprocess.run([sys.executable, s], cwd=BASE, stdout=fh, stderr=subprocess.STDOUT)
        status = "ok" if r.returncode == 0 else "FAILED"
        if r.returncode != 0 and s in LOCAL:
            status = "skipped (local data not available or failed; see log)"
        man["steps"].append(
            dict(
                script=s,
                status=status,
                seconds=round(time.time() - t0, 1),
                log=str(logp.relative_to(BASE)),
            )
        )
        print(f"{s}: {status} ({time.time() - t0:.0f} s)", flush=True)
        save()
        if r.returncode != 0 and s not in LOCAL:
            ok = False
            break
    if ok:
        try:
            from _r1_common import load_all

            df = load_all()
            man["sample"] = dict(
                N=len(df), groups={k: int(v) for k, v in df["k2"].value_counts().items()}
            )
        except Exception as e:
            man["sample"] = dict(error=str(e))
        outs = {}
        for p in sorted(AN.glob("S_R1*/*")):
            if p.is_file() and p.suffix in (".csv", ".md", ".png", ".docx"):
                outs[str(p.relative_to(AN))] = sha(p)
        man["outputs"] = outs
    man["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    man["all_ok"] = ok
    save()
    line = (
        f"| {run_id} | {man['pipeline'].get('pipeline_version')} {man['pipeline'].get('git_hash')}（{man['pipeline'].get('ica_algorithm')}） | {man['inputs']['psd_n']} | "
        f"{man.get('sample', {}).get('N', '–')} | {'全部成功' if ok else '失敗：' + man['steps'][-1]['script']} | work/runs/{run_id}/ |\n"
    )
    rl = BASE / "work" / "RUNLOG.md"
    if not rl.exists():
        rl.write_text(
            "# R1 分析的執行紀錄\n\n| 執行編號 | 前處理版本（演算法） | 頻譜檔數 | 樣本 N | 結果 | 紀錄位置 |\n|---|---|---|---|---|---|\n",
            encoding="utf-8",
        )
    with open(rl, "a", encoding="utf-8") as f:
        f.write(line)
    print("all done" if ok else "STOPPED", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
