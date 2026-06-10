function r = compute_data_rank(EEG, n_interpolated)
% COMPUTE_DATA_RANK  計算 EEG 資料之有效 rank,給 ICA PCA pre-reduction 用
%
% USAGE
%   r = compute_data_rank(EEG, n_interpolated)
%
% 邏輯(MNE 文檔標準):
%   rank = nbchan - 1 (avg ref 後) - n_interpolated
%
% 若 EEG 尚未做 avg ref(step1 / step2a 階段),則不扣 1。
% 此 helper 假設「將要做 ICA + 之後會 avg ref」之 step2,所以扣 1。
%
% 用法:
%   pop_runica(EEG, 'icatype', 'picard', 'pca', r, ...)
%
% v2.2 ─ 2026/5/27

if nargin < 2; n_interpolated = 0; end

% 基本 rank
r = EEG.nbchan;

% 扣 avg reference 將造成的 rank deficiency
% 注意:即使 step2 ICA 前還沒 avg ref,我們是要為 step2b 之後的 avg ref 預留 rank
r = r - 1;

% 扣每個被插補的通道(每個 interp 通道是 rank-deficient 之 線性組合)
r = r - n_interpolated;

% 安全下限
if r < 1
    r = 1;
end
end
