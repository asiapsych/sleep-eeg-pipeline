# 靜息態 EEG 處理與分析

[English](README.md) | **中文**

第 1–3 階段將 NeuroGuide 腦波轉成 EDF，完成前處理後再擷取頻譜特徵。這裡只提供程式碼，受試者資料由研究團隊保管。

*Sleep Medicine* 稿件使用的前處理，以及 ISI 題項的主成分分析與分群程式，見
[`04_manuscript_sleepmed_r1/`](04_manuscript_sleepmed_r1/README.md)。

---

## 處理步驟

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
├── CHANGELOG.md
├── 01_ng2edf/                # Python，將 NeuroGuide .ng 轉成 EDF
│   ├── ng2edf.py
│   ├── NG_FORMAT.md              # .ng 格式的位元組層級說明（英文）
│   ├── tests/                    # 以合成 .ng 檔進行的單元測試
│   └── README.md
├── 02_preprocess_eeglab/     # MATLAB，前處理與匯出 PSD
│   ├── run_pipeline_auto.m       # 全自動入口
│   ├── run_pipeline_manual.m     # 人工確認入口（論文採用）
│   ├── master_pipeline.m         # 流程協調（step1 … step5）
│   ├── step1 … step5、review_*_smart、config/、helpers/、chan_locs/
│   └── README.md
├── 03_features/              # Python，從 PSD 計算頻段功率與 FOOOF
└── 04_manuscript_sleepmed_r1/   # Sleep Medicine 稿件專用的程式
    ├── matlab_pipeline_v2.3/     # 分析資料實際使用的前處理版本
    ├── analysis/                 # ISI 題項的主成分分析與分群（Python）
    └── README.md
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
# 1) 將 NeuroGuide 檔案轉成 EDF   （加 --events 可另存每個檔案的錄製接縫表）
python 01_ng2edf/ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf
```

階段 1 補充：

- 轉檔程式依檔頭中記錄的區段位置（section offsets）與通道代碼表，以結構化方式解析 `.ng`，只輸出代碼 0–18 的 19 個 EEG 通道。
- 平均通道或檔案大小對不上時，該檔會報錯並略過，不會寫出錯位的資料。確定要轉時可加 `--lenient`。
- 檔尾每個取樣點各有一個事件字組。目前只見過 `0x00020000`，代表錄製中斷後接續的「接縫」（splice）：該點與下一點之間的訊號不連續，最可能是擷取時按了暫停再繼續。轉檔時，有接縫的檔案會各印一行警告，列出接縫時間；`--events` 另寫出 `<ID>_events.tsv`（欄位 sample、time_s、word_hex、kind、note）。
- EDF 仍是一般連續的 EDF，不標記接縫。分析時請把接縫當成段落邊界：epoch、頻譜視窗、濾波與偵測器都不可跨越接縫，接縫前後各捨棄至少 1 秒。
- 使用者手動放置的 Event Marker 存在哪裡、每個取樣點那一個位元組的意義，目前仍不清楚。
- 格式細節與尚未解開的問題見 [`01_ng2edf/NG_FORMAT.md`](01_ng2edf/NG_FORMAT.md)（英文）。
- 測試：`python -m unittest discover -s 01_ng2edf/tests -v`。

```matlab
% 2) MATLAB 前處理（互動式 MATLAB 中）
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
run_pipeline_manual('/path/to/raw_edf')   % 人工確認（論文模式）
% 或 run_pipeline_auto('/path/to/raw_edf') % 全自動
```

```bash
# 3) 從 PSD 擷取特徵
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

第 1 階段從輸入資料夾名稱取得去識別化代碼，寫入 EDF 的受試者欄位。轉檔程式不會複製 NeuroGuide 檔頭裡的身分資訊。

- 檔頭的受試者資訊區可能含姓名、出生日期、備註等，轉檔程式一律不解碼、不顯示、不寫出。
- 年齡欄位（舊版程式誤稱為 version）同樣不寫入任何輸出。
- EDF 的開始日期與時間取自 NeuroGuide 的施測日期與時間。

---

## 授權（License）

本 repo 為從研究室內部 EEG 分析 pipeline **去識別化／通用化**後整理出的研究重現程式碼，以 **BSD 3-Clause License** 釋出（見 `LICENSE`）。

> This public release was prepared by extracting, de-identifying, and generalizing the reusable analysis components from the laboratory's internal study-specific pipeline.

---

## 程式來源、去識別化與 AI 協助

研究流程由作者開發。公開前使用 Claude（Anthropic）拆出可重用的程式，並移除受試者編號、個別檔名及本機路徑。第 01–03 階段使用 Claude Opus 4.8，`04_manuscript_sleepmed_r1/` 使用 Claude Opus 5.5。編號和路徑由本機設定檔讀取，或以範例值替代。作者在首次發布前核對了去識別化後的程式。

Codex（OpenAI）協助整理公開版的 Python 主成分分析與分群程式。整理後重跑腳本，並與上一版比對數值輸出。

---

## 引用（Citation）

重複使用或修改這套腦波轉檔、前處理與頻譜分析程式時，請引用相關論文及本程式庫。

修訂稿使用的封存版本：v2.0.2，https://doi.org/10.5281/zenodo.24415381。

**建議的資料/程式碼可用性聲明（針對稿件）：**

> 用於腦電圖轉換、預處理、功率譜密度提取和特徵生成的去識別化和通用化分析程式碼可在 https://github.com/sleepaffect/insomnia-resting-eeg 取得。原始腦電圖記錄、參與者層面的臨床數據、資料庫文件和識別碼映射均不公開，因為本研究的倫理審批和參與者知情同意書並未包含參與者層面數據的公開共享，且這些數據受機構隱私限制。

---

## 維護者

**睡眠與情緒科學研究室**（Sleep and Affective Science Lab）
亞洲大學心理學系 臨床心理學組
網站：https://sleepaffect.org

## 致謝

本流程建構於 [EEGLAB](https://sccn.ucsd.edu/eeglab/)（前處理／ICA）、[ICLabel](https://github.com/sccn/ICLabel) 與 `clean_rawdata` plugins，以及 [FOOOF / specparam](https://fooof-tools.github.io/fooof/)（頻譜參數化）。
