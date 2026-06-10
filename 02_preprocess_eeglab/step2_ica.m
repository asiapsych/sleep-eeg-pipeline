function step2_ica(sid_list, opts)
% STEP2_ICA  CleanLine 60Hz + Rank-reduced ICA + ICLabel(僅標記)
%
% USAGE
%   step2_ica();
%   step2_ica({'sub01'});
%   step2_ica([], struct('overwrite', true));
%
% 流程(對齊 v2.2 §4):
%   1. 載入 step1.set
%   2. CleanLine 60Hz(取代 Notch,不破壞鄰近頻率)
%   3. 計算 data rank(扣 avg ref + n interp)
%   4. rng(set_subject_rng) → pop_runica(picard 或 runica,PCA=rank)
%   5. iclabel
%   6. 存 ica.set + 寫 iclabel.json + manifest + log
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
pipe  = config_pipeline();

fprintf('============================================================\n');
fprintf('  Step 2: CleanLine + Rank-reduced ICA + ICLabel\n');
fprintf('  Pipeline %s | %d subjects | ICA: %s\n', ...
        pipe.version, numel(sid_list), thr.ica.algorithm);
fprintf('============================================================\n\n');

t_start = tic;
n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};

for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.LAB_DATA, [sid '_step1.set']);
    out_path = fullfile(paths.LAB_DATA, [sid '_ica.set']);

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(in_path, 'file')
        fprintf('SKIP (no step1)\n');
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
        EEG = pop_loadset('filename', [sid '_step1.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);

        % 2. CleanLine 60Hz
        try
            EEG = pop_cleanline(EEG, ...
                'bandwidth',         thr.cleanline.bandwidth, ...
                'chanlist',          1:EEG.nbchan, ...
                'computepower',      0, ...
                'linefreqs',         thr.cleanline.line_freq, ...
                'normSpectrum',      0, ...
                'p',                 0.01, ...
                'pad',               2, ...
                'plotfigures',       0, ...
                'scanforlines',      thr.cleanline.scanForLines, ...
                'sigtype',           'Channels', ...
                'tau',               thr.cleanline.tau, ...
                'verb',              0, ...
                'winsize',           4, ...
                'winstep',           1);
        catch ME_cl
            % Fallback to Notch
            fprintf('(CleanLine failed: %s; fall back to Notch) ', ME_cl.message);
            EEG = pop_eegfiltnew(EEG, ...
                thr.cleanline.line_freq - 1, ...
                thr.cleanline.line_freq + 1, [], 1);
        end
        EEG = eeg_checkset(EEG);

        % 3. Data rank(估 n_interpolated 從 manifest)
        mf_path = fullfile(paths.MANIFEST_DIR, sprintf('%s.json', sid));
        n_interp = 0;
        if exist(mf_path, 'file')
            try
                mf = jsondecode(fileread(mf_path));
                if isfield(mf, 'step1') && isfield(mf.step1, 'n_bad')
                    n_interp = mf.step1.n_bad;
                end
            catch
            end
        end
        data_rank = compute_data_rank(EEG, n_interp);

        % 4. Set per-subject seed + ICA(根據演算法選 options)
        seed = set_subject_rng(sid);
        try
            switch lower(thr.ica.algorithm)
                case 'picard'
                    % Picard:用 EEGLAB 預設(maxiter=500, mode='standard')+ PCA
                    EEG = pop_runica(EEG, ...
                        'icatype', 'picard', ...
                        'maxiter', 500, ...
                        'mode', 'standard', ...
                        'pca', data_rank, ...
                        'interupt', 'off');
                case 'runica'
                    EEG = pop_runica(EEG, ...
                        'icatype', 'runica', ...
                        'extended', thr.ica.extended, ...
                        'pca', data_rank, ...
                        'interupt', 'off');
                otherwise
                    EEG = pop_runica(EEG, ...
                        'icatype', thr.ica.algorithm, ...
                        'pca', data_rank, ...
                        'interupt', 'off');
            end
        catch ME_ica
            if strcmp(thr.ica.algorithm, 'picard')
                fprintf('(picard failed: %s; fall back to runica) ', ME_ica.message);
                EEG = pop_runica(EEG, ...
                    'icatype', 'runica', ...
                    'extended', 1, ...
                    'pca', data_rank, ...
                    'interupt', 'off');
            else
                rethrow(ME_ica);
            end
        end
        EEG = eeg_checkset(EEG);
        n_comp = size(EEG.icaweights, 1);

        % 5. ICLabel
        EEG = iclabel(EEG);
        scores  = EEG.etc.ic_classification.ICLabel.classifications;
        classes = EEG.etc.ic_classification.ICLabel.classes;
        eye_col    = find(strcmpi(classes, 'Eye'));
        muscle_col = find(strcmpi(classes, 'Muscle'));

        % 自動 flag(尚不刪)
        auto_flag = [];
        flag_reason = {};
        for ci = 1:n_comp
            eye_s = scores(ci, eye_col);
            mus_s = scores(ci, muscle_col);
            if eye_s >= thr.iclabel.eye_thresh
                auto_flag(end+1) = ci; %#ok<AGROW>
                flag_reason{end+1} = sprintf('IC%d:Eye=%.2f', ci, eye_s);
            elseif mus_s >= thr.iclabel.muscle_thresh
                auto_flag(end+1) = ci; %#ok<AGROW>
                flag_reason{end+1} = sprintf('IC%d:Muscle=%.2f', ci, mus_s);
            end
        end
        EEG.etc.iclabel_auto_flag = auto_flag;
        EEG.etc.iclabel_v22       = struct('seed', seed, 'rank', data_rank);

        % 6. Save
        EEG = pop_saveset(EEG, 'filename', [sid '_ica.set'], 'filepath', paths.LAB_DATA);

        % 7. iclabel.json
        ic_json = struct();
        ic_json.subject       = sid;
        ic_json.n_components  = n_comp;
        ic_json.classes       = {classes};
        ic_json.scores        = scores;
        ic_json.auto_flag     = auto_flag;
        ic_json.flag_reason   = {flag_reason};
        ic_json.seed          = seed;
        ic_json.rank          = data_rank;
        ic_json.ica_algorithm = thr.ica.algorithm;
        fid = fopen(fullfile(paths.DECISIONS, [sid '_iclabel.json']), 'w');
        fprintf(fid, '%s', jsonencode(ic_json, 'PrettyPrint', true));
        fclose(fid);

        % 8. Manifest
        info = struct();
        info.input        = in_path;
        info.output       = out_path;
        info.n_components = n_comp;
        info.data_rank    = data_rank;
        info.n_interp     = n_interp;
        info.seed         = seed;
        info.ica_algorithm = thr.ica.algorithm;
        info.cleanline_freq = thr.cleanline.line_freq;
        info.iclabel_auto_flag = auto_flag;
        info.iclabel_flag_count = numel(auto_flag);
        manifest_writer(sid, 'step2', info);

        % 9. Log
        log_writer('step2', sid, ...
            {'n_components','data_rank','n_interp','seed','ica_algo','n_flagged','flag_reason','status','note'}, ...
            {n_comp, data_rank, n_interp, seed, thr.ica.algorithm, ...
             numel(auto_flag), strjoin(flag_reason, ';'), 'ok', ''});

        fprintf('OK (n_comp=%d, rank=%d, flagged=%d)\n', n_comp, data_rank, numel(auto_flag));
        n_ok = n_ok + 1;
    catch ME
        log_writer('step2', sid, ...
            {'n_components','data_rank','n_interp','seed','ica_algo','n_flagged','flag_reason','status','note'}, ...
            {NaN, NaN, NaN, NaN, '', 0, '', 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('  Step 2 done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Elapsed: %.1f min\n', elapsed/60);
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('  Next step: review_ica_smart() then step2b_apply_ica()\n');
fprintf('============================================================\n');
end
