function [bad_labels, scores] = detect_bad_channels(EEG)
% DETECT_BAD_CHANNELS  Multi-metric bad channel detection(v2.2 簡化版)
%
% USAGE
%   [bad, scores] = detect_bad_channels(EEG);
%
% 4 個指標,任 ≥2 觸發 = bad:
%   1. clean_rawdata ChannelCriterion < 0.6
%   2. Flatline > 5s
%   3. Line noise (60Hz) SNR > 10 dB
%   4. Channel kurtosis robust z > 5
%
% RETURN
%   bad_labels : cell array of channel labels
%   scores     : struct;每個 channel 之 4 個指標
%
% v2.2 ─ 2026/5/27

thr = config_thresholds();

n_ch = EEG.nbchan;
labels = {EEG.chanlocs.labels};

% 初始化 score table
scores = struct();
scores.labels = labels;
scores.fail_corr = false(1, n_ch);
scores.fail_flat = false(1, n_ch);
scores.fail_line = false(1, n_ch);
scores.fail_kurt = false(1, n_ch);

% ── 指標 1+2:clean_rawdata(corr + flatline)──
% clean_rawdata(EEG, FlatlineSec, HighpassMinHz, ChanCorr, LineNoise, BurstCrit, WindowCrit)
EEG_test = clean_rawdata(EEG, thr.badch.flatline_sec, -1, ...
                          thr.badch.corr_thresh, -1, -1, -1);
remaining = {EEG_test.chanlocs.labels};
removed_by_clean = setdiff(labels, remaining);
for k = 1:length(removed_by_clean)
    idx = find(strcmp(labels, removed_by_clean{k}));
    if ~isempty(idx)
        % clean_rawdata 不細分 corr vs flat,兩個都算 fail_corr
        scores.fail_corr(idx) = true;
        scores.fail_flat(idx) = true;
    end
end

% ── 指標 3:Line noise SNR ──
data = EEG.data;
if ndims(data) == 3
    data = reshape(data, n_ch, []);
end
for ci = 1:n_ch
    [pxx, f] = pwelch(data(ci, :), hamming(min(round(EEG.srate*2), size(data,2))), ...
                       [], [], EEG.srate);
    % 60Hz 附近 (58-62) 與相鄰段 (50-55 / 65-70) 比較
    p_60 = mean(pxx((f >= 58) & (f <= 62)));
    p_bg = mean(pxx(((f >= 50) & (f <= 55)) | ((f >= 65) & (f <= 70))));
    if p_bg > 0
        snr_db = 10 * log10(p_60 / p_bg);
        if snr_db > thr.badch.line_noise_snr_db
            scores.fail_line(ci) = true;
        end
    end
end

% ── 指標 4:Channel kurtosis robust z ──
ch_kurt = zeros(1, n_ch);
for ci = 1:n_ch
    ch_kurt(ci) = kurtosis(data(ci, :));
end
% Robust z 用 median + MAD
med = median(ch_kurt);
mad_val = median(abs(ch_kurt - med));
if mad_val > 0
    rz = (ch_kurt - med) / (1.4826 * mad_val);
else
    rz = zeros(size(ch_kurt));
end
scores.fail_kurt = abs(rz) > thr.badch.kurtosis_robust_z;
scores.kurtosis_z = rz;

% ── 投票:任 ≥2 觸發 = bad ──
vote_count = double(scores.fail_corr) + double(scores.fail_flat) + ...
             double(scores.fail_line) + double(scores.fail_kurt);
% 注意:fail_corr + fail_flat 是同源(clean_rawdata),所以 corr+flat 算 1 票
% 重新計算:clean_rawdata 1 票 + line 1 票 + kurt 1 票 = max 3
vote_count2 = double(scores.fail_corr | scores.fail_flat) + ...
              double(scores.fail_line) + double(scores.fail_kurt);
scores.vote_count = vote_count2;

is_bad = vote_count2 >= 2;
bad_labels = labels(is_bad);

scores.is_bad = is_bad;
end
