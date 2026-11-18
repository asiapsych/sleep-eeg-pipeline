# 01_ng2edf — NeuroGuide (.ng) → EDF

> 語言：**繁體中文** · [English](README.md)

針對 NeuroGuide 專有 `.ng` 二進位格式，以逆向工程開發的讀取器，
並具備符合標準的 EDF 寫入器。NeuroGuide 未提供公開的轉換工具或
規格。本讀取器遵循的位元組配置記載於
**[NG_FORMAT.md](NG_FORMAT.md)**。

檔頭依結構解析：

- 依據檔案中儲存的區段位移，定位 float32 資料與尾端
  區段。
- 區段 B 提供各資料區塊的取樣率、取樣點數與電極
  代碼。
- 僅保留 EEG 區塊（電極代碼 0..18）。
- 一致性檢查會停止轉換配置不符的檔案，
  避免寫入錯位的資料。

## 環境需求

- Python 3.8+
- `numpy`（請參閱最上層的 `requirements.txt`）

## 使用方式

```bash
python ng2edf.py --input /path/to/raw_ng --output /path/to/raw_edf

# Preview without writing
python ng2edf.py --input ./data/raw_ng --output ./data/raw_edf --dry-run

# Only process case IDs with a given prefix
python ng2edf.py -i ./data/raw_ng -o ./data/raw_edf --accept-prefix S

# Also write <ID>_events.tsv with the events (recording splices) of each file
python ng2edf.py -i ./data/raw_ng -o ./data/raw_edf --events
```

| 選項 | 作用 |
|---|---|
| `--events` | 在每個 EDF 旁另行寫入 `<ID>_events.tsv`（請參閱[事件：錄製接縫](#events-recording-splices)）。 |
| `--lenient` | 對未通過檔案大小或平均通道檢查的檔案，仍進行轉換並發出警告，而非回報錯誤。 |

若所有檔案皆轉換成功，結束代碼為 0；若有任何檔案失敗，則為 2。失敗的
檔案會列為 `ERROR` 行，批次處理仍會繼續。

## 輸入結構

`--input` 下的每個直接子資料夾代表一名受試者，並存放一個 `.ng`：

```
raw_ng/
├── S001=/        →  S001.edf
│   └── xxx.ng
├── S002=/        →  S002.edf
│   └── yyy.ng
└── ...
```

個案 ID 為資料夾名稱開頭符合 `[A-Za-z]+[0-9]+` 的字串，並轉為大寫
（`S001=` → `S001`、`A012X` → `A012`）。若資料夾內有多個 `.ng` 檔案，則使用
排序後的第一個檔案。

## 通道

EEG 通道依 NeuroGuide 電極代碼（0..18）選取，並依
19 通道的 10-20 順序寫入：

```
FP1 FP2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 Fz Cz Pz
```

檔案中的兩個虛擬通道不會寫入 EDF：

- 連結耳參考通道（代碼 55），其值全為零；
- 平均通道（代碼 56），僅用於完整性檢查。

導程配置的參考方式（例如 LinkEars）會以後綴形式附加至通道名稱
（`FP1-LE`）。下游 MATLAB 階段會進一步標準化通道（移除
後綴，轉換為現代 10-20 命名 T7/T8/P7/P8）。

NeuroGuide 儲存為 NaN 的通道（例如已移除的電極）會
以 0 寫入，並發出警告。請在下游將其標記為不良通道。

## 檢查

每個檔案皆依 [NG_FORMAT.md](NG_FORMAT.md) §10 的配置進行檢查。

- **一律視為錯誤**：
  - 缺少魔術值；
  - 區段位移不一致；
  - 區段 A 或 B 未在下一個區段位移處結束；
  - 檔案在 EEG 資料內遭截斷。
- **視為錯誤，或在使用 `--lenient` 時發出警告**：
  - 檔案大小不等於尾端位移 + 5 × 取樣點數；
  - 平均通道（代碼 56）與 EEG 通道的平均值
    不符。這表示通道表遭到誤讀。
- **僅發出警告**：未知的電極代碼、未知的導程配置名稱、
  非預期的前導常數、非零填補值、尾端各取樣點的非零位元組
  （連同其數量）、無法辨識的事件字組，以及未位於
  6 個取樣點區塊末端的接縫。
- **每個含接縫的檔案皆會發出警告**：列出接縫時間，例如
  `recording discontinuity at t=12.343 s; do not let analysis windows span it`。

## 事件：錄製接縫

`.ng` 檔案的尾端區段為每個取樣點儲存一個 `u32` 事件字組。
目前觀察到唯一的非零字組為 `0x00020000`。它標示一處**接縫**：錄製資料
在該取樣點與下一個取樣點之間不連續，最可能是因為
資料擷取曾暫停後再繼續。相關證據記載於
[NG_FORMAT.md](NG_FORMAT.md) §8；已在 384 筆錄製資料上驗證。接縫
並非使用者放置的事件標記。NeuroGuide 如何儲存這些標記，目前仍
不清楚，因為檢查過的檔案中未出現這類標記。

轉換工具會針對每個含接縫的檔案印出一則警告，並列出
接縫時間。`--events` 另會將每個非零事件字組各寫成一列至
`<ID>_events.tsv`。即使錄製資料沒有事件，仍會寫入此檔案，
此時檔案僅含標頭列。

| 欄位 | 意義 |
|---|---|
| `sample` | 事件字組的取樣點索引，從 0 起算 |
| `time_s` | `sample / sample rate`，自第一個取樣點起算的秒數 |
| `word_hex` | 事件字組，例如 `0x00020000` |
| `kind` | `splice` 用於 `0x00020000`，`unknown` 用於其他值 |
| `note` | 簡短說明，例如 `recording discontinuous after this sample` |

**分析注意事項。** 輸出的 EDF 為一般的連續 EDF，不會標記
接縫。請將每個接縫視為區段邊界：不連續處位於
`time_s` 與 `time_s + 1/fs` 之間。

- 勿讓分段、頻譜視窗、濾波器或事件偵測器跨越接縫。
- 請分別處理接縫兩側，並捨棄接縫前後各至少 1 s
  的資料。接縫後會出現短暫的暫態訊號。
- 若接縫源於暫停，接縫後的時間也會比實際時鐘時間
  落後一段暫停的時間，而檔案並未儲存該時長。

早期此格式讀取器所提供的標記位置，約比正確位置晚了
整筆錄製時長的四分之一。NG_FORMAT.md §11 提供換算方式。

## 隱私

- 輸出 EDF 的*病人識別*欄位設為
  **去識別化的個案 ID**（取自資料夾名稱），**不會**設為任何
  可能儲存在 NeuroGuide 檔頭中的識別碼。請確認輸入
  資料夾名稱已完成去識別化。
- NeuroGuide 的受試者資訊區段可能包含姓名、出生日期
  及備註。轉換工具絕不解碼、印出或寫入這些字串。
- 轉換工具會讀取年齡欄位，但絕不寫入。EDF 的開始
  日期／時間取自 NeuroGuide 的檢測日期／時間。
- 請將 EDF 檔案與事件表皆視為參與者資料。兩者皆由
  儲存庫的 `.gitignore` 排除（`*.edf`、`*_events.tsv`）。

## 測試

測試僅使用合成的 `.ng` 檔案（由 `tests/ngsynth.py` 在
暫存資料夾中建立），其中包含將 EDF 輸出與
前一版轉換工具逐位元組比對的檢查。請從儲存庫根目錄執行：

```bash
python -m unittest discover -s 01_ng2edf/tests -v
```

## 備註

- 物理值最小值／最大值對所有通道使用單一全域範圍（單一增益）；
  進行對絕對振幅敏感的分析時，請留意此點。
- 此格式透過逆向工程解析，並已針對一種
  擷取配置（DSI-24、19 通道、300 Hz、連結耳參考）的錄製資料完成驗證。若使用其他
  配置，請檢查警告，並使用 EDF
  檢視器或 MNE 抽查數個轉換後的檔案。
