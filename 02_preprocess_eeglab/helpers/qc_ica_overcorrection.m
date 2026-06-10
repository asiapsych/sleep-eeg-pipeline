function result = qc_ica_overcorrection(EEG_pre, EEG_post, sid)
% QC_ICA_OVERCORRECTION  自動偵測 ICA 是否過度校正(2024 best practice)
%
% USAGE
%   r = qc_ica_overcorrection(EEG_step1, EEG_step2, 'sub01');
%
% 邏輯(Bailey 2024 / Coelli 2024 / Delorme 2023):
%   - 計算每個 primary band 之 PSD %change
%   - 任一頻段(尤其 Beta 13-30)|% change| > 20% → 紅旗
%   - Brain 訊號通常不應被 ICA 大幅移除,過大 change = over-correction
%
% RETURN
%   result.pct_change_per_band  : 各頻段 % change
%   result.beta_pct_change      : Beta 頻段 % change(主指標)
%   result.flag                 : true 若 Beta |% change| > thr
%   result.warn_msg             : 描述
%
% v2.2 ─ 2026/5/27

bands = config_bands();
thr   = config_thresholds();

% 對兩個 EEG 算 BRAIN-avg PSD(連續訊號)
[psd_pre,  f_pre]  = compute_brain_avg_psd(EEG_pre);
[psd_post, f_post] = compute_brain_avg_psd(EEG_post);

primary_names = fieldnames(bands.primary);
n_band = numel(primary_names);

pct_change = zeros(1, n_band);
for i = 1:n_band
    rng = bands.primary.(primary_names{i});
    mp = (f_pre >= rng(1)) & (f_pre <= rng(2));
    mq = (f_post >= rng(1)) & (f_post <= rng(2));
    pb = trapz(f_pre(mp), psd_pre(mp));
    pa = trapz(f_post(mq), psd_post(mq));
    if pb > 0
        pct_change(i) = (pa - pb) / pb * 100;
    end
end

% Beta band 主指標
beta_idx = find(strcmp(primary_names, 'Beta'));
beta_change = pct_change(beta_idx);

flag = abs(beta_change) > thr.overcorr.warn_pct_change;

result = struct();
result.subject = sid;
result.band_names = primary_names;
result.pct_change_per_band = pct_change;
result.beta_pct_change = beta_change;
result.flag = flag;
result.threshold = thr.overcorr.warn_pct_change;
if flag
    result.warn_msg = sprintf( ...
        '[OVER-CORRECTION FLAG] %s: Beta 改變 %.1f%% (門檻 ±%d%%) ─ ICA 可能誤刪 Beta', ...
        sid, beta_change, thr.overcorr.warn_pct_change);
else
    result.warn_msg = '';
end
end


function [psd, f] = compute_brain_avg_psd(EEG)
thr = config_thresholds();
data = EEG.data;
if ndims(data) == 3
    data = reshape(data, EEG.nbchan, []);
end
nperseg = min(thr.psd.win_samples, size(data, 2));
psd_all = [];
for ci = 1:EEG.nbchan
    [pxx, f] = pwelch(data(ci, :), hamming(nperseg), 0, nperseg, EEG.srate);
    if isempty(psd_all)
        psd_all = zeros(EEG.nbchan, length(pxx));
    end
    psd_all(ci, :) = pxx';
end
psd = mean(psd_all, 1);
end
