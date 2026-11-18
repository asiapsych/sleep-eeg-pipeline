# 02_preprocess_eeglab — EDF → 前處理 → PSD

> 語言：**繁體中文** · [English](README.md)

以 EEGLAB 為基礎、全程以腳本執行的靜息態 EEG 前處理流程，從 EDF
產生可供分析的 PSD（以及選用的頻帶功率）。以函式為基礎，具備集中式
設定、確定性種子、版本註記、品質管制與日誌記錄。

## 環境需求

- MATLAB R2020a 或更新版本
- [EEGLAB](https://sccn.ucsd.edu/eeglab/)（2021+），已安裝並加入 MATLAB 路徑
- EEGLAB 外掛：`clean_rawdata`、`ICLabel`、`firfilt`；`Picard` 為選用
  （無法使用時改用 `runica`）

> **未**隨附 EEGLAB。請先執行 `addpath('/path/to/eeglab'); eeglab nogui;`；
> 流程會自動偵測 EEGLAB 的位置。

## 兩種模式：自動化與人工確認

首先在 MATLAB 載入 EEGLAB：

```matlab
addpath('/path/to/eeglab');  eeglab nogui;
cd 02_preprocess_eeglab
```

### 模式 A — 全自動（`run_pipeline_auto`）

```matlab
run_pipeline_auto('/path/to/raw_edf');     % step1→step5, no human input
```

所有 ICA 成分均由 ICLabel 判定，分段則依規則剔除；無需
人工介入。最適合快速處理大量資料。

### 模式 B — 人工確認（`run_pipeline_manual`，**論文採用此模式**）

```matlab
run_pipeline_manual('/path/to/raw_edf');   % machine runs heavy steps; human confirms
```

電腦執行前處理與 ICA，接著**由人工確認兩項有歧義的
判定** — 此為論文採用的工作流程：

| 確認項目 | 腳本 | 確認內容 |
|--------------------|--------|------------------|
| ICA 成分 | `review_ica_smart` | 僅確認**中等信心度**的成分（ICLabel 最高分數位於審查區間內）；高信心度者自動移除，低信心度者自動保留。顯示地形圖／頻譜／活化訊號；按鍵 `r`=移除 `k`=保留 `s`=略過受試者 `q`=儲存並離開 |
| 分段 | `review_epochs_smart` | 審查低品質受試者的分段剔除結果 |

> 兩種模式採用完全相同的 step1–step5 邏輯與參數；差別僅在
> ICA 成分與分段由自動確認或人工確認。
> 人工判定會寫入 `outputs/decisions/<ID>_ica_decisions.json`，並在
> 重新執行時沿用。手動模式必須在**互動式** MATLAB 工作階段中執行。

### 進階：分階段／自訂執行

```matlab
% Manual mode with explicit start/stop steps
opts = struct('overwrite', false, 'skip_review', false, ...
              'start_from', 'step1', 'stop_after', 'step5');
master_pipeline([], opts);

% Re-run only some steps for one subject
master_pipeline({'sub01'}, struct('start_from','step3','stop_after','step4'));
```

## 步驟

| 步驟 | 檔案 | 執行內容 |
|------|------|--------------|
| 1 | `step1_preprocess.m` | 載入 EDF、標準化通道（移除 -LE、T3/4/5/6→T7/8/P7/P8）、套用 10-20 位置、0.5–70 Hz 帶通濾波、多指標不良通道偵測（相關性／平坦訊號／電源雜訊／峰度，≥2 項觸發）、球面插值 |
| 2 | `step2_ica.m` | CleanLine 60 Hz（無法使用時改用陷波濾波）、降秩 ICA（Picard，無法使用時改用 runica）、ICLabel 分類 |
| 2b | `step2b_apply_ica.m` | 套用移除判定（Eye ≥ 0.70 / Muscle ≥ 0.80）、平均參考、ICA 過度校正品質管制 |
| 3 | `step3_epoch_reject.m` | 3 s 分段、三層剔除（autorej／振幅 150 µV／跳變 100 µV）、抽取 60 個分段 |
| 4 | `step4_export_psd.m` | Welch PSD（Hamming，無重疊），跨分段取平均 → `psd_<ID>.mat` |
| 5 | `step5_export_bandpower.m` | （選用）由 PSD 計算絕對／相對頻帶功率 → CSV |

## 輸入／輸出

- **輸入**：`raw_edf/*.edf`（受試者 ID = EDF 基本檔名）。使用
  `../01_ng2edf` 將 NeuroGuide 資料轉換為這些檔案。
- **輸出**（皆位於 `outputs/` 下，由 git 排除）：
  - `psd_mat/psd_<ID>.mat` — 供 `../03_features` 使用
  - `band_power/` — 選用的 MATLAB 端頻帶功率
  - `qc_pdf/`、`qc_png/`、`logs/`、`manifests/`、`methods/`

> 建議：將 `../03_features`（Python）視為傳統頻帶功率的唯一權威來源，
> 與 FOOOF 共用相同的 PSD 與頻帶定義。
> `step5` 的 MATLAB 頻帶功率為選用／供交叉核對。

## 設定

| 檔案 | 內容 |
|------|----------|
| `config/config_paths.m` | 路徑（CHAN_LOC 指向隨附的 `chan_locs/`；EEGLAB_DIR 自動偵測） |
| `config/config_thresholds.m` | 所有閾值（濾波器、不良通道、ICA、ICLabel、分段、PSD、A/B/C/F 評等） |
| `config/config_bands.m` | 頻帶定義 |
| `config/config_pipeline.m` | 版本、19 通道、排除清單、隨機種子 |
| `config/config_user.m` | **各專案的覆寫設定（由 git 排除）**。複製 `config_user.m.example` 後編輯 |

## 其他工具

- `master_qc_report.m` — 每名受試者的 5 頁品質管制 PDF（A/B/C/F 評等）
- `master_methods_export.m` — 自動產生 COBIDAS 風格的方法段落草稿
- `population_qc_summary.m` — 全體受試者的品質管制概覽
- `install_picard.m` — 安裝 Picard 外掛（請先載入 EEGLAB）
- `review_ica_smart.m` / `review_epochs_smart.m` — 半自動人工審查

## 可再現性

- 各受試者的確定性種子（受試者 ID + 全域種子的雜湊值）
- `version_stamp` 記錄 git 雜湊值 + EEGLAB／外掛版本 + 時間戳記
- 每名受試者的 JSON 資訊清單 + 每個步驟的 CSV 日誌
- 注意：`runica` 並非完全確定性；Picard 較穩定，但仍
  取決於 BLAS 的實作
