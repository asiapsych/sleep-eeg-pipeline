function step6_export_psd_laplacian(sid_list)
% STEP6_EXPORT_PSD_LAPLACIAN  表面拉普拉斯轉換後的頻譜（v2.3 新增；肌電的敏感度分析用）
%
% USAGE
%   step6_export_psd_laplacian();            % 所有有 60 個 epoch 的受試者
%   step6_export_psd_laplacian({'S000001'});
%
% 讀 <SID>_step3.set（60 個 3 秒 epoch，與 step4 用的是同一份資料），
% 先對每個時間點做球面樣條表面拉普拉斯（helpers/surface_laplacian_perrin.m；
% stiffness 4、lambda 1e-5、Legendre 50 項），再照 step4 完全相同的設定算 Welch 頻譜
% （Hamming、win/nfft 取自 config_thresholds、無重疊、跨 epoch 平均）。
% 存成 outputs/psd_mat_laplacian/psd_<SID>.mat，欄位與 step4 相同：
%   freqs, psd [n_ch × n_freq], ch_names, n_ep, sfreq, subject
% 另存 T（轉換矩陣）與 xyz（電極座標）供核對。
%
% 主分析不使用這一步的輸出。
%
% v2.3 ─ 2026/11/17

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
paths = config_paths(); thr = config_thresholds();
out_dir = fullfile(paths.PIPELINE_ROOT, 'outputs', 'psd_mat_laplacian');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

n_ok = 0; n_skip = 0; n_fail = 0;
for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.LAB_DATA, [sid '_step3.set']);
    out_path = fullfile(out_dir, sprintf('psd_%s.mat', sid));
    if ~exist(in_path, 'file') || exist(out_path, 'file'); n_skip = n_skip + 1; continue; end
    try
        EEG = pop_loadset('filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);
        if EEG.trials ~= thr.epoch.target_count; n_skip = n_skip + 1; continue; end
        sfreq = EEG.srate; n_ch = EEG.nbchan; n_ep = EEG.trials; ch_names = {EEG.chanlocs.labels};

        xyz = [[EEG.chanlocs.X]', [EEG.chanlocs.Y]', [EEG.chanlocs.Z]'];
        if size(xyz, 1) ~= n_ch || any(isnan(xyz(:)))
            error('%s: 電極座標不完整，無法計算表面拉普拉斯', sid);
        end
        T = surface_laplacian_perrin(xyz);
        lap = reshape(T * reshape(double(EEG.data), n_ch, []), n_ch, EEG.pnts, n_ep);

        win = hamming(thr.psd.win_samples);
        [~, f] = pwelch(squeeze(lap(1, :, 1)), win, thr.psd.overlap, thr.psd.nfft, sfreq);
        f_mask = (f >= thr.psd.freq_min) & (f <= thr.psd.freq_max); freqs = f(f_mask);
        psd = zeros(n_ch, numel(freqs));
        for ch = 1:n_ch
            acc = zeros(1, numel(freqs));
            for tr = 1:n_ep
                pxx = pwelch(squeeze(lap(ch, :, tr)), win, thr.psd.overlap, thr.psd.nfft, sfreq);
                acc = acc + pxx(f_mask)';
            end
            psd(ch, :) = acc / n_ep;
        end

        subject = sid; %#ok<NASGU>
        save(out_path, 'freqs', 'psd', 'ch_names', 'n_ep', 'sfreq', 'subject', 'T', 'xyz', '-v7');
        log_writer('step6', sid, {'n_ep','n_ch','sfreq','n_freqs','status','note'}, ...
                   {n_ep, n_ch, sfreq, numel(freqs), 'ok', 'surface Laplacian (Perrin; m=4, lambda=1e-5, 50 terms)'});
        n_ok = n_ok + 1;
        if mod(n_ok, 25) == 0; fprintf('  step6: %d done\n', n_ok); end
    catch ME
        log_writer('step6', sid, {'n_ep','n_ch','sfreq','n_freqs','status','note'}, ...
                   {NaN, NaN, NaN, NaN, 'error', ME.message});
        fprintf('  step6 FAIL %s: %s\n', sid, ME.message);
        n_fail = n_fail + 1;
    end
end
fprintf('Step 6 done: %d ok, %d skip, %d fail / total %d\n  Output: %s\n', n_ok, n_skip, n_fail, numel(sid_list), out_dir);
end
