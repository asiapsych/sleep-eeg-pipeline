function step1_preprocess(sid_list, opts)
% STEP1_PREPROCESS  Load EDF + Bandpass + Multi-metric Bad Channel + Interp
%
% USAGE
%   step1_preprocess();                    % 全部 subjects(default)
%   step1_preprocess({'sub01','sub02'});
%   step1_preprocess([], struct('overwrite', true));
%
% 流程(對齊 v2.2 設計筆記 §4):
%   1. 載入 EDF + 19-channel 標準化(load_edf_clean.m)
%   2. Bandpass 0.5-70 Hz(pop_eegfiltnew)
%   3. Multi-metric bad channel detection(detect_bad_channels.m)
%   4. Spherical spline interpolation(pop_interp)
%   5. QC PSD plot(qc_psd_plot.m)
%   6. 存 step1.set + manifest + log
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
pipe  = config_pipeline();

fprintf('============================================================\n');
fprintf('  Step 1: Load + Bandpass + Bad Channel + Interp\n');
fprintf('  Pipeline %s | %d subjects | overwrite=%d\n', ...
        pipe.version, numel(sid_list), opts.overwrite);
fprintf('============================================================\n\n');

t_start = tic;
n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};

for i = 1:numel(sid_list)
    sid = sid_list{i};
    edf_path = fullfile(paths.RAW_EDF_DIR, [sid '.edf']);
    out_path = fullfile(paths.LAB_DATA, [sid '_step1.set']);

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(edf_path, 'file')
        fprintf('SKIP (no EDF)\n');
        n_skip = n_skip + 1;
        continue;
    end
    if exist(out_path, 'file') && ~opts.overwrite
        fprintf('SKIP (already done)\n');
        n_skip = n_skip + 1;
        continue;
    end

    try
        % 1. Load + clean channels
        EEG = load_edf_clean(edf_path);
        n_chan_in = EEG.nbchan;
        srate     = EEG.srate;
        duration  = EEG.xmax;
        EEG_raw = EEG;  % 保留 raw 供 QC

        % 2. Bandpass
        EEG = pop_eegfiltnew(EEG, thr.bandpass.highpass, thr.bandpass.lowpass);
        EEG = eeg_checkset(EEG);

        % 3. Multi-metric bad channel detection
        [bad_labels, bad_scores] = detect_bad_channels(EEG);
        n_bad = numel(bad_labels);
        bad_str = strjoin(bad_labels, ';');

        % 4. Interpolation(若有 bad)
        if n_bad > 0
            % 先記原 chanlocs
            orig_chanlocs = EEG.chanlocs;
            bad_idx = find(ismember({EEG.chanlocs.labels}, bad_labels));
            EEG = pop_select(EEG, 'nochannel', bad_idx);
            EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
            EEG = eeg_checkset(EEG);
        end

        % 5. Save
        EEG.setname = [sid '_step1'];
        EEG = pop_saveset(EEG, 'filename', [sid '_step1.set'], 'filepath', paths.LAB_DATA);

        % 6. QC plot(bandpass 前後 + bad ch topography)
        try
            qc_psd_plot(EEG_raw, EEG, sid, 'step1', paths.QC_PNG);
        catch ME_qc
            fprintf('(qc plot fail: %s) ', ME_qc.message);
        end

        % 7. Manifest
        info = struct();
        info.input              = edf_path;
        info.output             = out_path;
        info.n_chan_in          = n_chan_in;
        info.srate              = srate;
        info.duration_s         = duration;
        info.bandpass           = [thr.bandpass.highpass, thr.bandpass.lowpass];
        info.bad_channels       = bad_labels;
        info.n_bad              = n_bad;
        info.bad_vote_count     = bad_scores.vote_count;
        info.bad_kurt_z         = bad_scores.kurtosis_z;
        info.bad_fail_corr      = bad_scores.fail_corr;
        info.bad_fail_flat      = bad_scores.fail_flat;
        info.bad_fail_line      = bad_scores.fail_line;
        info.bad_fail_kurt      = bad_scores.fail_kurt;
        manifest_writer(sid, 'step1', info);

        % 8. Log
        log_writer('step1', sid, ...
            {'n_chan_in','srate','duration_s','bad_channels','n_bad','status','note'}, ...
            {n_chan_in, srate, duration, bad_str, n_bad, 'ok', ''});

        fprintf('OK (bad=%d:%s)\n', n_bad, bad_str);
        n_ok = n_ok + 1;
    catch ME
        log_writer('step1', sid, ...
            {'n_chan_in','srate','duration_s','bad_channels','n_bad','status','note'}, ...
            {NaN, NaN, NaN, '', 0, 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('  Step 1 done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Elapsed: %.1f min\n', elapsed/60);
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('  Output: %s\n', paths.LAB_DATA);
fprintf('============================================================\n');
end
