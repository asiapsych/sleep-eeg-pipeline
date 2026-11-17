function step5_export_epoch_psd(sid_list)
% STEP5_EXPORT_EPOCH_PSD  匯出逐 epoch 的頻譜（v2.3 新增；清醒度的逐段分析用）
%
% USAGE
%   step5_export_epoch_psd();            % 所有有 60 個 epoch 的受試者
%   step5_export_epoch_psd({'S000001'});
%
% 讀 <SID>_step3.set（60 個 3 秒 epoch），對每個 epoch、每個電極各算一次 Welch 頻譜
% （設定與 step4 完全相同：Hamming、win/nfft 取自 config_thresholds、無重疊），
% 不做平均，存成 outputs/psd_epochs/psdep_<SID>.mat：
%   psd_ep   [n_ep × n_ch × n_freq]
%   freqs, ch_names, n_ep, sfreq, epoch_index（step3 抽樣後的順序）
% 對 60 個 epoch 取平均即等於 step4 的 psd。
%
% v2.3 ─ 2026/11/16

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
paths = config_paths(); thr = config_thresholds();
out_dir = fullfile(paths.PIPELINE_ROOT, 'outputs', 'psd_epochs');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

n_ok = 0; n_skip = 0;
for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.LAB_DATA, [sid '_step3.set']);
    out_path = fullfile(out_dir, sprintf('psdep_%s.mat', sid));
    if ~exist(in_path, 'file') || exist(out_path, 'file'); n_skip = n_skip + 1; continue; end
    try
        EEG = pop_loadset('filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);
        if EEG.trials ~= thr.epoch.target_count; n_skip = n_skip + 1; continue; end
        sfreq = EEG.srate; n_ch = EEG.nbchan; n_ep = EEG.trials; ch_names = {EEG.chanlocs.labels};
        win = hamming(thr.psd.win_samples);
        [~, f] = pwelch(squeeze(EEG.data(1, :, 1)), win, thr.psd.overlap, thr.psd.nfft, sfreq);
        f_mask = (f >= thr.psd.freq_min) & (f <= thr.psd.freq_max); freqs = f(f_mask);
        psd_ep = zeros(n_ep, n_ch, numel(freqs));
        for ch = 1:n_ch
            for tr = 1:n_ep
                pxx = pwelch(squeeze(EEG.data(ch, :, tr)), win, thr.psd.overlap, thr.psd.nfft, sfreq);
                psd_ep(tr, ch, :) = pxx(f_mask);
            end
        end
        epoch_index = 1:n_ep; subject = sid; %#ok<NASGU>
        save(out_path, 'freqs', 'psd_ep', 'ch_names', 'n_ep', 'sfreq', 'subject', 'epoch_index', '-v7');
        log_writer('step5', sid, {'n_ep','n_ch','n_freqs','status','note'}, {n_ep, n_ch, numel(freqs), 'ok', ''});
        n_ok = n_ok + 1;
    catch ME
        fprintf('FAIL %s: %s\n', sid, ME.message);
        log_writer('step5', sid, {'n_ep','n_ch','n_freqs','status','note'}, {NaN, NaN, NaN, 'fail', ME.message});
    end
end
fprintf('step5 epoch PSD: %d exported, %d skipped\n', n_ok, n_skip);
end
