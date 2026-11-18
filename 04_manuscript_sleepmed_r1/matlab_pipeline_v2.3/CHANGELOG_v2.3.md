# v2.3（2026-11-16）

v2.2 的最終執行（2026-05-28）裡，獨立成分分析與電源雜訊移除兩步的呼叫都失敗，程式靜靜改用備援方法
（extended Infomax、59–61 Hz 帶拒濾波），而紀錄檔寫的仍是設定值（picard）。v2.3 修掉這個問題並重新執行。

## 更動

| 檔案 | 更動 | 原因 |
|---|---|---|
| `step2_ica.m` | 呼叫 PICARD 時拿掉 `'interupt'` 參數 | picard 不認得這個參數，v2.2 因此每一位都失敗 |
| `step2_ica.m` | 呼叫 CleanLine 前把外掛自帶的輔助函式加進路徑 | v2.2 找不到 `hlp_varargin2struct`，每一位都失敗 |
| `step2_ica.m` | 拿掉兩處備援（PICARD → runica、CleanLine → notch）；失敗就記為該受試者失敗 | 不再靜靜換方法 |
| `step1_preprocess.m`、`step2b_apply_ica.m` | 批次執行時不畫每人的品管圖 | MATLAB R2026b 無視窗模式下繪圖輔助程序會當掉並卡住主程序 |
| `config/config_pipeline.m` | 版本改為 v2.3 | |
| `step5_export_epoch_psd.m` | 新增：匯出逐 epoch 的頻譜（設定與 step4 相同，不取平均） | 清醒度的逐段分析需要 |
| `run_all_v23.m` | 新的啟動腳本：從原始 EDF 跑到頻譜匯出，全程寫紀錄 | |
| `step6_export_psd_laplacian.m`、`helpers/surface_laplacian_perrin.m`、`run_step6_v23.m` | 新增（2026-11-17）：同一批 60 個 epoch 先做球面樣條表面拉普拉斯（Perrin et al., 1989；stiffness 4、lambda 1e-5、Legendre 50 項），再用與 step4 相同的設定算頻譜，存到 `outputs/psd_mat_laplacian/`。轉換矩陣以模擬資料和 MNE-Python 的 `compute_current_source_density` 比對過，逐電極相關為 1，只差固定倍數（頭半徑平方）| 肌電的敏感度分析需要；主分析不使用這一步的輸出 |

## 沒有更動的

濾波（0.5–70 Hz）、壞通道判準、ICLabel 門檻（眼動 ≥ .70、肌肉 ≥ .80）、平均參考、3 秒 epoch、三層 epoch 剔除、
每人隨機抽 60 個 epoch（每人固定亂數種子）、Welch 設定。

## 與 v2.2 輸出的關係

v2.2 的輸出（extended Infomax＋notch）保留在 `matlab_pipeline_v2/outputs/`，作為對照。
兩者的逐人一致性與主要結果的比較見分析程式庫的 `task_r1p_picard_comparison.py`。
