function step3_epoch_reject(sid_list, opts)
% STEP3_EPOCH_REJECT  3s epoch + 3-layer reject + sample 60 (per-subject seed)
%
% USAGE
%   step3_epoch_reject();
%   step3_epoch_reject({'sub01'});
%   step3_epoch_reject([], struct('overwrite', true));
%
% 流程:
%   1. 載入 step2.set
%   2. 切 3s 固定 epoch(eeg_regepochs)
%   3. Layer A: pop_autorej(振幅 + kurtosis 機率)
%   4. Layer B: 振幅硬門檻 150uV
%   5. Layer C: jump 100uV
%   6. 採樣 60 個(per-subject seed);若 n_good < 60 全取 + 加 review queue
%   7. 存 step3.set + manifest + log + epoch_quality.json
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
pipe  = config_pipeline();

fprintf('============================================================\n');
fprintf('  Step 3: Epoch + 3-Layer Reject + Sample %d\n', thr.epoch.target_count);
fprintf('  Pipeline %s | %d subjects\n', pipe.version, numel(sid_list));
fprintf('============================================================\n\n');

t_start = tic;
n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};
need_review_list = {};   % epoch <60 或 reject% > 25% 之個案

for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.LAB_DATA, [sid '_step2.set']);
    out_path = fullfile(paths.LAB_DATA, [sid '_step3.set']);

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(in_path, 'file')
        fprintf('SKIP (no step2)\n');
        n_skip = n_skip + 1;
        continue;
    end
    if exist(out_path, 'file') && ~opts.overwrite
        fprintf('SKIP (done)\n');
        n_skip = n_skip + 1;
        continue;
    end

    try
        % 1. Load
        EEG = pop_loadset('filename', [sid '_step2.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);

        % 2. Epoch
        EEG = eeg_regepochs(EEG, 'recurrence', thr.epoch.length_sec, ...
            'limits', [0 thr.epoch.length_sec], 'rmbase', NaN);
        EEG = eeg_checkset(EEG);
        n_total = EEG.trials;
        if n_total == 0
            error('No epoch could be cut');
        end

        % 3. Layer A: autorej
        [EEG, ~] = pop_autorej(EEG, 'nogui', 'on', ...
            'threshold',  thr.epoch.autorej_threshold, ...
            'startprob',  thr.epoch.autorej_startprob, ...
            'maxrej',     thr.epoch.autorej_maxrej);
        EEG = eeg_checkset(EEG);
        n_autorej = EEG.trials;

        % 4. Layer B: amp hard threshold
        data_uv = EEG.data;
        amp_bad = squeeze(any(any(abs(data_uv) > thr.epoch.amp_thresh_uv, 1), 2));
        EEG = pop_select(EEG, 'trial', find(~amp_bad));
        EEG = eeg_checkset(EEG);
        n_amp = EEG.trials;

        % 5. Layer C: jump
        diffs = diff(EEG.data, 1, 2);
        jump_bad = squeeze(any(any(abs(diffs) > thr.epoch.jump_thresh_uv, 1), 2));
        EEG = pop_select(EEG, 'trial', find(~jump_bad));
        EEG = eeg_checkset(EEG);
        n_jump = EEG.trials;

        if n_jump == 0
            error('All epochs rejected');
        end

        % 6. Sample 60 with per-subject seed
        seed = set_subject_rng(sid);
        target = thr.epoch.target_count;
        if n_jump < target
            % 全取 + 標 review
            sel_idx = 1:n_jump;
            below_target = true;
            need_review_list{end+1} = sprintf('%s(%d epochs)', sid, n_jump); %#ok<AGROW>
        else
            sel_idx = sort(randperm(n_jump, target));
            below_target = false;
        end
        EEG = pop_select(EEG, 'trial', sel_idx);
        EEG = eeg_checkset(EEG);
        n_final = EEG.trials;

        % Artifact reject rate should be computed before the fixed-count sampling.
        % n_final mixes true artifact rejection with downsampling to 60 epochs.
        reject_rate = 1 - n_jump / n_total;
        sampling_drop_rate = 1 - n_final / n_total;
        if reject_rate > thr.epoch.review_max_reject
            need_review_list{end+1} = sprintf('%s(artifact reject%%=%.0f)', sid, reject_rate*100); %#ok<AGROW>
        end

        % 7. Save
        EEG.setname = [sid '_step3'];
        EEG = pop_saveset(EEG, 'filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);

        % 8. epoch_quality.json
        eq = struct();
        eq.subject       = sid;
        eq.n_total       = n_total;
        eq.n_after_layerA = n_autorej;
        eq.n_after_layerB = n_amp;
        eq.n_after_layerC = n_jump;
        eq.n_final       = n_final;
        eq.reject_rate   = reject_rate;
        eq.sampling_drop_rate = sampling_drop_rate;
        eq.below_target  = below_target;
        eq.seed          = seed;
        eq.selected_idx  = sel_idx;
        fid = fopen(fullfile(paths.DECISIONS, [sid '_epoch_quality.json']), 'w');
        fprintf(fid, '%s', jsonencode(eq, 'PrettyPrint', true));
        fclose(fid);

        % 9. Manifest
        info = struct();
        info.input             = in_path;
        info.output            = out_path;
        info.epoch_length_sec  = thr.epoch.length_sec;
        info.n_total           = n_total;
        info.n_after_layerA    = n_autorej;
        info.n_after_layerB    = n_amp;
        info.n_after_layerC    = n_jump;
        info.n_final           = n_final;
        info.reject_rate       = reject_rate;
        info.sampling_drop_rate = sampling_drop_rate;
        info.below_target      = below_target;
        info.seed              = seed;
        info.target_count      = target;
        manifest_writer(sid, 'step3', info);

        % 10. Log
        log_writer('step3', sid, ...
            {'n_total','n_layerA','n_layerB','n_layerC','n_final','reject_rate','below_target','seed','status','note'}, ...
            {n_total, n_autorej, n_amp, n_jump, n_final, reject_rate, below_target, seed, 'ok', ''});

        fprintf('OK (epochs %d→A%d→B%d→C%d → %d, artifact reject=%.0f%%)\n', ...
                n_total, n_autorej, n_amp, n_jump, n_final, reject_rate*100);
        n_ok = n_ok + 1;
    catch ME
        log_writer('step3', sid, ...
            {'n_total','n_layerA','n_layerB','n_layerC','n_final','reject_rate','below_target','seed','status','note'}, ...
            {NaN, NaN, NaN, NaN, NaN, NaN, false, NaN, 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('  Step 3 done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Elapsed: %.1f min\n', elapsed/60);
fprintf('  Need review (low epoch or high reject%%): %d\n', numel(need_review_list));
if ~isempty(need_review_list)
    review_path = fullfile(paths.LOGS, 'step3_need_review.txt');
    fid = fopen(review_path, 'w');
    fprintf(fid, '# Subjects needing review_epochs_smart\n');
    for k = 1:numel(need_review_list)
        fprintf(fid, '%s\n', need_review_list{k});
    end
    fclose(fid);
    fprintf('  Written: %s\n', review_path);
end
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('============================================================\n');
end
