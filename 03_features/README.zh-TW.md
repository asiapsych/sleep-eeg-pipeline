# 03_features — PSD → 傳統功率 + FOOOF（週期性／非週期性）

> 語言：**繁體中文** · [English](README.md)

從階段 02 匯出的 `psd_<ID>.mat` 檔案，使用
**相同的 PSD 與相同的頻帶定義**計算兩類特徵：

1. **傳統頻帶功率** — 絕對／相對功率 + 常用比值（TBR、ABR、…）
2. **週期性／非週期性分解（FOOOF / specparam）** — 非週期性指數
   與偏移量，以及各頻帶的週期性峰值（中心頻率／功率／頻寬）

此腳本**僅進行特徵擷取**。分組、臨床測量
與統計檢定屬於各研究的非公開分析，不包含於此。

## 環境需求

- Python 3.8+
- `numpy / scipy / pandas / fooof`（請參閱最上層的 `requirements.txt`；
  FOOOF 1.x 需要 numpy < 2）

## 使用方式

```bash
python compute_features.py --psd-dir /path/to/outputs/psd_mat --out-dir /path/to/features

# Skip specific subjects; customize the fit frequency range
python compute_features.py --psd-dir ./psd_mat --out-dir ./features \
    --exclude sub01 sub02 --freq-min 1 --freq-max 45
```

## 輸入

`psd_<ID>.mat`（階段 02 `step4_export_psd` 的輸出），每個檔案包含：

| 變數 | 形狀 | 意義 |
|----------|-------|---------|
| `freqs` | `[n_freq]` | 頻率軸（Hz） |
| `psd` | `[n_ch × n_freq]` | 功率頻譜密度 |
| `ch_names` | `[n_ch]` | 通道標籤 |
| `n_ep` | 純量（選用） | 取平均的分段數 |

## 輸出（CSV）

| 檔案 | 內容 |
|------|----------|
| `traditional_power_subject.csv` | 每名受試者 × 各頻帶的絕對／相對功率（全腦 + 5 個腦區）+ 比值 |
| `specparam_channel.csv` | 每名受試者 × 各通道：指數／偏移量／峰值／R² |
| `specparam_subject.csv` | specparam 特徵的全腦／各腦區平均值 |
| `fit_quality.csv` | 每名受試者的 FOOOF 擬合 R²／誤差 |

## 頻帶（文獻標準）

- **主要頻帶**（連續、不重疊；相對功率的分母）：
  Delta 1–4、Theta 4–8、Alpha 8–13、Beta 13–30、Gamma 30–45 Hz
- **衍生子頻帶**（以占總量的比例呈現，不列入分母）：
  Alpha1 8–10、Alpha2 10–13、Beta1 13–18、Beta2 18–24、Beta3 24–30、HighBeta 30–40 Hz
- 相對功率 = 頻帶絕對功率 / 5 個主要頻帶的總和
  （= 單一連續的 1–45 Hz 積分）× 100
- 比值（TBR、ABR、…）使用絕對功率（總量會相消）

## 腦區（19 通道，現代 10-20）

| 腦區 | 通道 |
|--------|----------|
| 額區 | FP1 FP2 F3 F4 F7 F8 FZ |
| 中央區 | C3 C4 CZ |
| 頂區 | P3 P4 PZ |
| 顳區 | T7 T8 P7 P8 |
| 枕區 | O1 O2 |

> 若輸入 PSD 使用舊式顳區命名（T3/T4/T5/T6），載入時會將其標準化為
> 現代命名（T7/T8/P7/P8），因此不影響腦區計算。

## FOOOF 設定

`freq_range=(1,45)`、`aperiodic_mode="fixed"`、`peak_width_limits=(1,12)`、
`max_n_peaks=8`、`min_peak_height=0.1`、`peak_threshold=2.0`。
