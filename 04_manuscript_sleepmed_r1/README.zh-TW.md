# Sleep Medicine 稿件（修訂版）的程式

> 語言：**繁體中文** · [English](README.md)

本資料夾為失眠症狀嚴重度與清醒週期性 beta 活動之修訂稿（稿件編號 SLEEP-D-26-02563）專用的程式。本倉庫的 `01_ng2edf`、`02_preprocess_eeglab` 與 `03_features` 是通用、可重複使用的流程；本資料夾收錄稿件實際執行的前處理，以及失眠嚴重度量表（Insomnia Severity Index, ISI）題項的主成分分析與分群程式。

**不含任何受試者資料。** 依研究倫理審查的規定，原始記錄、問卷資料、受試者編號與所有分析輸出均由研究團隊保存。因此，沒有資料無法重新執行這些程式，但每一個步驟都可以在此閱讀。

```
04_manuscript_sleepmed_r1/
├── matlab_pipeline_v2.3/   # MATLAB + EEGLAB：分析資料實際使用的前處理
└── analysis/               # Python：ISI 題項的主成分分析與分群
```

MATLAB 程式的註解使用繁體中文，Python 程式的註解使用英文。

## `matlab_pipeline_v2.3/`：前處理

自動執行，所有參與者使用同一組固定門檻（`config/config_thresholds.m`）：帶通濾波、壞通道偵測與內插、CleanLine、以 PICARD 估計 ICA、依 ICLabel 移除成分、平均參考、切分 epoch 並剔除，最後輸出頻譜。自動標記的成分與保留的 epoch 另由受過訓練的研究人員目視確認（`review_ica_smart.m`、`review_epochs_smart.m`）。

| 檔案 | 用途 |
|---|---|
| `run_all_v23.m` | 入口：對所有記錄執行步驟 1–5 |
| `run_step6_v23.m` | 入口：只執行步驟 6 |
| `step1_preprocess.m` … `step3_epoch_reject.m` | 濾波、通道、ICA、epoch |
| `step4_export_psd.m` | 60 個 epoch 平均的功率譜 |
| `step5_export_epoch_psd.m` | 每個 epoch 的功率譜 |
| `step6_export_psd_laplacian.m` | 球面樣條表面拉普拉斯轉換後的頻譜 |
| `validation_runica_vs_picard.m`、`validation_test_retest.m` | ICA 演算法與重複執行一致性的檢查 |
| `config/` | 路徑、流程設定、門檻、頻段 |
| `CHANGELOG_v2.3.md` | 相對於前一版流程的變更 |

`config/config_paths.m` 為範本：資料位置由環境變數 `EEG_RAW_EDF_DIR`、`EEG_RAW_NG_DIR` 與 `EEGLAB_DIR` 讀取。前處理之前即排除的參與者列於本機檔案 `config/excluded_subjects.txt`，不包含在倉庫中。

## `analysis/`：主成分分析與分群

所有程序僅使用 ISI 總分 9 分以上參與者的七個 ISI 題項。各腳本逐支執行。

| 腳本 | 內容 |
|---|---|
| `task_r1w_pc_interpretation.py` | 七個 ISI 題項的主成分分析：負荷量、解釋變異、平行分析（2,000 組置換資料）與負荷量的重抽樣一致性（2,000 次） |
| `task_e2_bootstrap_selectK.py` | 二至六群的預測回原樣本重抽樣穩定度（1,000 次） |
| `task_r1b_k234_comparison.py` | 二至六群的內部效度指標（silhouette、Calinski–Harabasz、Davies–Bouldin），並列 `task_e2` 的重抽樣穩定度 |
| `task_r1e_clustering_robustness.py` | 兩群解在九種分群方法（K-means、四種高斯混合模型、三種階層式方法、潛在類別分析）之間的一致性，以及 gap statistic |
| `_dynamic_n.py`、`_r1_common.py` | 共用程式：樣本定義、資料讀取，以及 K-means 兩組（低嚴重度、高嚴重度） |

`task_r1b_k234_comparison.py` 與 `task_r1e_clustering_robustness.py` 會讀取 `task_e2_bootstrap_selectK.py` 的重抽樣結果，請先執行後者。

### 輸入

- 問卷資料的 SQLite 資料庫（ISI 題項與總分、年齡、性別、教育年數）；讀取的欄位列於 `_r1_common.py` 與 `_dynamic_n.py`。
- 頻譜特徵表（週期性與非週期性參數、相對頻段功率；見第 03 階段），用來界定具可用 EEG 的參與者。

資料庫位置與排除參與者的編號由只留在本機的 `local_settings.json` 讀取；格式見 `local_settings.example.json`。若設定環境變數 `PSYML_DB`，以環境變數為準。

### 環境

Python 3.12；套件版本見 `analysis/requirements.txt`。
