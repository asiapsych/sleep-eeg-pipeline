function step4_export_psd(sid_list, opts)
% STEP4_EXPORT_PSD  PSD .mat export per subject(對齊 task_f / Bai 2025)
%
% USAGE
%   step4_export_psd();
%   step4_export_psd({'S000001'});
%   step4_export_psd([], struct('overwrite', true));
%
% 計算:
%   每 channel 每 epoch 跑 pwelch(Hamming 3s window, no overlap, nfft=900)
%   → 跨 epoch 平均
%   freq range output: 0.5-70 Hz, resolution 0.333 Hz
%
% 輸出:psd_S<id>.mat
%   fields: freqs (1×F), psd (n_ch × F), ch_names (cell), n_ep, sfreq, subject
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
pipe  = config_pipeline();

fprintf('============================================================\n');
fprintf('  Step 4: PSD Export(.mat)\n');
fprintf('  Pipeline %s | %d subjects\n', pipe.version, numel(sid_list));
fprintf('  pwelch: %d-sample Hamming, no overlap, nfft=%d, range %.1f-%.0f Hz\n', ...
        thr.psd.win_samples, thr.psd.nfft, thr.psd.freq_min, thr.psd.freq_max);
fprintf('============================================================\n\n');

t_start = tic;
n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};

for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.LAB_DATA, [sid '_step3.set']);
    out_path = fullfile(paths.PSD_OUT, sprintf('psd_%s.mat', sid));

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(in_path, 'file')
        fprintf('SKIP (no step3)\n');
        n_skip = n_skip + 1;
        continue;
    end
    if exist(out_path, 'file') && ~opts.overwrite
        fprintf('SKIP (done)\n');
        n_skip = n_skip + 1;
        continue;
    end

    try
        EEG = pop_loadset('filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);
        sfreq    = EEG.srate;
        n_ch     = EEG.nbchan;
        n_ep     = EEG.trials;
        ch_names = {EEG.chanlocs.labels};

        % Strict analysis rule: only exactly 60 three-second epochs are exported.
        % Low-epoch step3 files may remain for audit, but must not enter Python.
        if n_ep ~= thr.epoch.target_count
            error('Strict n=%d rule: %s has n_ep=%d; PSD export skipped', ...
                  thr.epoch.target_count, sid, n_ep);
        end

        % Pre-compute freq mask
        win = hamming(thr.psd.win_samples);
        [pxx_tmp, f] = pwelch(squeeze(EEG.data(1, :, 1)), win, thr.psd.overlap, ...
                              thr.psd.nfft, sfreq);
        f_mask = (f >= thr.psd.freq_min) & (f <= thr.psd.freq_max);
        freqs = f(f_mask);
        n_f = length(freqs);

        psd = zeros(n_ch, n_f);
        for ch = 1:n_ch
            psd_ch = zeros(n_ep, n_f);
            for tr = 1:n_ep
                x = squeeze(EEG.data(ch, :, tr));
                [pxx, ~] = pwelch(x, win, thr.psd.overlap, thr.psd.nfft, sfreq);
                psd_ch(tr, :) = pxx(f_mask)';
            end
            psd(ch, :) = mean(psd_ch, 1);
        end

        subject = sid; %#ok<NASGU>
        save(out_path, 'freqs', 'psd', 'ch_names', 'n_ep', 'sfreq', 'subject', '-v7');

        % Manifest
        info = struct();
        info.input    = in_path;
        info.output   = out_path;
        info.n_ep     = n_ep;
        info.n_ch     = n_ch;
        info.sfreq    = sfreq;
        info.n_freqs  = n_f;
        info.win_samples = thr.psd.win_samples;
        info.nfft     = thr.psd.nfft;
        manifest_writer(sid, 'step4', info);

        % Log
        log_writer('step4', sid, ...
            {'n_ep','n_ch','sfreq','n_freqs','status','note'}, ...
            {n_ep, n_ch, sfreq, n_f, 'ok', ''});

        fprintf('OK (n_ep=%d, n_f=%d)\n', n_ep, n_f);
        n_ok = n_ok + 1;
    catch ME
        log_writer('step4', sid, ...
            {'n_ep','n_ch','sfreq','n_freqs','status','note'}, ...
            {NaN, NaN, NaN, NaN, 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('  Step 4 done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Elapsed: %.1f min\n', elapsed/60);
fprintf('  Output: %s\n', paths.PSD_OUT);
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('============================================================\n');
end
