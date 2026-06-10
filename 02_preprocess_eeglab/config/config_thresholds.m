function thr = config_thresholds()
% CONFIG_THRESHOLDS  Centralized thresholds for the generic pipeline.
%
% USAGE
%   t = config_thresholds();
%   t.bandpass.highpass        % 0.5
%   t.iclabel.eye_thresh       % 0.70
%
thr = struct();

% ── Bandpass ──
thr.bandpass.highpass = 0.5;
thr.bandpass.lowpass  = 70.0;

% ── Notch / CleanLine ──
thr.cleanline.line_freq   = 60;
thr.cleanline.bandwidth   = 2;
thr.cleanline.scanForLines= true;
thr.cleanline.tau         = 100;

% ── Bad channel detection(multi-metric;任 ≥2 觸發 = bad)──
thr.badch.flatline_sec       = 5;       % flatline > 5s
thr.badch.corr_thresh        = 0.60;    % clean_rawdata corr
thr.badch.line_noise_snr_db  = 10;      % 60Hz peak SNR(dB) > 10
thr.badch.kurtosis_robust_z  = 5;       % channel kurtosis robust z

% ── ICA ──
thr.ica.algorithm    = 'picard';   % 'picard'(SOTA,已裝 v2.0)或 'runica'(備援)
thr.ica.extended     = 1;          % runica 用(picard 不接此參數)
thr.ica.picard_mode  = 'standard'; % picard 用('standard' 或 'ortho')
thr.ica.picard_maxiter = 500;
thr.ica.global_seed  = 42;         % 結合 per-subject hash 作 seed

% ── ICLabel(自動移除規則)──
thr.iclabel.eye_thresh        = 0.70;
thr.iclabel.muscle_thresh     = 0.80;
% Smart review 範圍(這區間需人工 review;範圍外自動 confirm)
thr.iclabel.review_min        = 0.40;
thr.iclabel.review_max        = 0.85;

% ── Over-correction QC(A7)──
thr.overcorr.beta_band             = [13 30];
thr.overcorr.warn_pct_change       = 20;    % Beta band % change > 20% → 紅旗
thr.overcorr.flag_subject_pct      = 5;     % 全人群中紅旗超過此比例 → pipeline-level warn

% ── Epoch ──
thr.epoch.length_sec       = 3.0;
thr.epoch.target_count     = 60;
thr.epoch.strict_target_count = false;  % false: export PSD even when fewer epochs survive.

% ── Three-layer reject ──
thr.epoch.autorej_threshold = 150;   % Layer A µV
thr.epoch.autorej_startprob = 5;     % Layer A probability SD
thr.epoch.autorej_maxrej    = 0.30;  % Layer A max 30%
thr.epoch.amp_thresh_uv     = 150;   % Layer B µV
thr.epoch.jump_thresh_uv    = 100;   % Layer C µV
thr.epoch.review_min_count  = 60;    % < N → enter review queue
thr.epoch.review_max_reject = 0.25;  % reject% > 25% → enter review queue

% ── PSD pwelch ──
thr.psd.win_sec      = 3.0;
thr.psd.win_samples  = [];           % []: auto round(win_sec * sampling rate)
thr.psd.nfft         = [];           % []: use win_samples
thr.psd.overlap      = 0;
thr.psd.freq_min     = 0.5;
thr.psd.freq_max     = 70;
thr.psd.compute_bandpower = true;

% ── Quality grading ──
thr.qc.grade_A_max_bad_ch     = 2;
thr.qc.grade_A_max_ic_removed = 7;
thr.qc.grade_A_min_epoch      = 60;
thr.qc.grade_B_max_bad_ch     = 4;
thr.qc.grade_B_max_ic_removed = 10;
thr.qc.grade_B_min_epoch      = 50;
% C = 之間;F = 任一指標超過 B 上限或無法完成 pipeline
end
