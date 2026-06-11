# 靜息態 EEG 處理 Pipeline（NeuroGuide → 特徵）

[English](README.md) | **中文**

從 NeuroGuide 原始腦波檔，一路處理到可供統計分析的**傳統頻段功率**與**週期性／非週期性（aperiodic / periodic, FOOOF）**特徵表的可重現 pipeline。

本 repo **只包含方法（程式碼）**，不含任何受試者資料。下游的分組、臨床量表與統計分析屬於各研究的私有部分，不在此開源範圍內。

---

## Pipeline 總覽

三段式線性流程，依序執行：

```
  原始 .ng (NeuroGuide)
        │
        ▼
┌─────────────────────────┐
│ 01_ng2edf  (Python)     │   逆向解析 NeuroGuide 二進位格式，輸出標準 EDF
└─────────────────────────┘
        │  *.edf
        ▼
┌─────────────────────────┐
│ 02_preprocess_eeglab    │   濾波 → 壞通道偵測/插補 → CleanLine → ICA(ICLabel)
│   (MATLAB + EEGLAB)     │   → 平均參考 → 切 epoch → 拒絕 → 匯出 PSD
└─────────────────────────┘
        │  psd_<ID>.mat
        ▼
┌─────────────────────────┐
│ 03_features  (Python)   │   傳統頻段功率（絕對/相對/比值）
│                         │   + FOOOF 非週期性 exponent/offset + 週期性峰值
└─────────────────────────┘
        │
        ▼
   特徵表 (CSV)  ──►  （下游統計分析，不在本 repo）
```

| 階段 | 語言 | 輸入 | 輸出 |
|------|------|------|------|
| `01_ng2edf` | Python | `.ng` 資料夾 | `.edf` |
| `02_preprocess_eeglab` | MATLAB + EEGLAB | `.edf` | `psd_<ID>.mat`、band power CSV |
| `03_features` | Python | `psd_*.mat` | 特徵表 CSV（傳統功率 + FOOOF） |

各階段詳見子資料夾的 `README.md`。

---

## 倉庫結構

```
.
├── README.md                 # 英文（預設）
├── README.zh-TW.md           # 本檔（中文）
├── LICENSE                   # BSD 3-Clause
├── requirements.txt          # 鎖版本的 Python 相依套件
├── .gitignore                # 排除所有資料與快取
├── 01_ng2edf/                # Python — NeuroGuide .ng → EDF
├── 02_preprocess_eeglab/     # MATLAB — EDF → 前處理 → PSD (+ band power)
│   ├── run_pipeline_auto.m       # 全自動入口
│   ├── run_pipeline_manual.m     # 人工確認入口（論文採用）
│   ├── master_pipeline.m         # 流程協調（step1 … step5）
│   ├── step1 … step5、review_*_smart、config/、helpers/、chan_locs/
│   └── README.md
└── 03_features/              # Python — PSD → 傳統功率 + FOOOF
```

---

## 需求

### Python（階段 1、3）

```bash
pip install -r requirements.txt
```

- 階段 1（轉檔）只需 `numpy`
- 階段 3（特徵）需 `numpy / scipy / pandas / fooof`（FOOOF 1.x 需 numpy < 2）

### MATLAB（階段 2）

- MATLAB R2020a 以上
- [EEGLAB](https://sccn.ucsd.edu/eeglab/)（2021+），需另外安裝並加入 path
- plugins：`clean_rawdata`、`ICLabel`、`firfilt`；`Picard` 可選（沒有時自動回退 `runica`）
- 本 repo **不內含 EEGLAB**；請先 `addpath('/path/to/eeglab'); eeglab nogui;`，pipeline 會自動偵測位置

---

## 快速開始

```bash
# 1) NeuroGuide → EDF
python 01_ng2edf/ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf
```

```matlab
% 2) MATLAB 前處理（互動式 MATLAB 中）
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
run_pipeline_manual('/path/to/raw_edf')   % 人工確認（論文模式）
% 或 run_pipeline_auto('/path/to/raw_edf') % 全自動
```

```bash
# 3) PSD → 特徵
python 03_features/compute_features.py \
    --psd-dir /path/to/outputs/psd_mat --out-dir /path/to/features
```

---

## 前處理：兩種模式（全自動 vs. 人工確認）

階段 2 提供兩個對等的一鍵入口，共用完全相同的 step1–step5 邏輯與參數，差別只在**模糊判斷**怎麼決定：

- **`run_pipeline_auto`** — 所有 ICA 成分與 epoch 全自動決定（ICLabel + 規則拒絕），適合大批快速處理。
- **`run_pipeline_manual`**（論文採用）— 機器跑完重步驟後，由人工確認兩處模糊判斷：中信心度的 ICA 成分（`review_ica_smart`，顯示 topoplot／頻譜／activation）與低品質個案的 epoch 拒絕（`review_epochs_smart`）。高/低信心自動處理；人工決定寫入 `outputs/decisions/<ID>_ica_decisions.json`，重跑沿用。

---

## 設計原則

- **可重現**：固定亂數種子；每位受試者用確定性 seed（hash subject ID + global seed）；版本戳記（git hash + 工具版本）；完整 per-subject JSON manifest 與 per-step CSV log。
- **單一功率來源**：傳統頻段功率與 FOOOF 都在階段 3 從**同一份 PSD、同一套頻段定義**計算，避免不同腳本頻段邊界不一致。
- **方法匯出**：階段 2 可執行 `master_methods_export` 自動產生 COBIDAS 規範的方法段落草稿。

### 頻段定義（文獻標準）

- **主頻段**（連續非重疊，進相對功率分母）：Delta 1–4、Theta 4–8、Alpha 8–13、Beta 13–30、Gamma 30–45 Hz
- **衍生子頻段**（只算占比，不進分母）：Alpha1 8–10、Alpha2 10–13、Beta1 13–18、Beta2 18–24、Beta3 24–30、HighBeta 30–40 Hz
- 相對功率 = 該頻段絕對功率 / 主頻段總和（= 1–45 Hz 連續積分）× 100

---

## 資料與隱私

本 repo **不包含任何**原始 EEG、受試者層級臨床資料、資料庫檔、ID 對照表或研究專屬 subject code。受試者層級資料受機構倫理審查、受試者同意與隱私法規限制，不公開。

Pipeline 在設計上即注重隱私：階段 1 轉檔時，輸出 EDF 的受試者欄位寫入**去識別化代碼**（依輸入資料夾名），不會帶入 NeuroGuide 檔頭中可能存在的原始識別資訊。

---

## 授權（License）

本 repo 為從研究室內部 EEG 分析 pipeline **去識別化／通用化**後整理出的研究重現程式碼，以 **BSD 3-Clause License** 釋出（見 `LICENSE`）。

> This public release was prepared by extracting, de-identifying, and generalizing the reusable analysis components from the laboratory's internal study-specific pipeline.

---

## 程式碼來源與 AI 協助

本公開倉庫是研究室內部、特定研究專屬分析 pipeline 的**去識別化與通用化改寫版**。原始程式碼包含受試者層級的檔案路徑、ID 對照表與資料庫耦合，無法公開分享。將可重用元件抽取、去識別化並通用化為本 release 的工作，由 **Claude（Claude Opus 4.8, Anthropic）** 協助完成。釋出之程式碼經作者審閱，並依 BSD 3-Clause License 以「現狀（as is）」提供，不附任何明示或默示之擔保（見 `LICENSE`）。

---

## 引用（Citation）

發布的程式碼支援腦電圖轉換、預處理、功率譜密度提取和特徵生成工作流程的可複現性。重複使用或修改程式碼時，請引用相關論文和此程式碼庫。

**建議的資料/程式碼可用性聲明（針對稿件）：**

> 用於腦電圖轉換、預處理、功率譜密度提取和特徵生成的去識別化和通用化分析程式碼可在 https://github.com/asiapsych/sleep-eeg-pipeline 取得。原始腦電圖記錄、參與者層面的臨床數據、資料庫文件和識別碼映射均不公開，因為本研究的倫理審批和參與者知情同意書並未包含參與者層面數據的公開共享，且這些數據受機構隱私限制。

---

## 維護者

**睡眠與情緒科學研究室**（Sleep and Affective Science Lab）
亞洲大學心理學系 臨床心理學組
網站：https://chang-lab-psychology-website-one.vercel.app/

## 致謝

本流程建構於 [EEGLAB](https://sccn.ucsd.edu/eeglab/)（前處理／ICA）、[ICLabel](https://github.com/sccn/ICLabel) 與 `clean_rawdata` plugins，以及 [FOOOF / specparam](https://fooof-tools.github.io/fooof/)（頻譜參數化）。
