function step2b_apply_ica(sid_list, opts)
% STEP2B_APPLY_ICA  Apply ICA decisions + Average Reference + Over-correction QC
%
% USAGE
%   step2b_apply_ica();
%   step2b_apply_ica({'sub01'});
%   step2b_apply_ica([], struct('overwrite', true));
%
% 流程:
%   1. 載入 ica.set + decisions.json
%   2. pop_subcomp(remove final_remove)
%   3. 平均參考 pop_reref(EEG, [])
%   4. 存 step2.set
%   5. Over-correction QC(step1 vs step2 PSD)─ A7 重點
%   6. Manifest + log + QC png
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
pipe  = config_pipeline();

fprintf('============================================================\n');
fprintf('  Step 2b: Apply ICA Decisions + Avg Ref + Over-correction QC\n');
fprintf('  Pipeline %s | %d subjects\n', pipe.version, numel(sid_list));
fprintf('============================================================\n\n');

t_start = tic;
n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};
overcorr_flags = {};   % 累積有 over-correction 紅旗者

for i = 1:numel(sid_list)
    sid = sid_list{i};
    ica_path  = fullfile(paths.LAB_DATA, [sid '_ica.set']);
    dec_path  = fullfile(paths.DECISIONS, [sid '_ica_decisions.json']);
    step1_path = fullfile(paths.LAB_DATA, [sid '_step1.set']);
    out_path  = fullfile(paths.LAB_DATA, [sid '_step2.set']);

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(ica_path, 'file')
        fprintf('SKIP (no ica.set)\n');
        n_skip = n_skip + 1;
        continue;
    end
    if ~exist(dec_path, 'file')
        fprintf('SKIP (no decision)\n');
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
        EEG = pop_loadset('filename', [sid '_ica.set'], 'filepath', paths.LAB_DATA);
        EEG = eeg_checkset(EEG);

        dec = jsondecode(fileread(dec_path));
        final_remove = dec.final_remove(:)';

        % 2. Remove IC
        if ~isempty(final_remove)
            EEG = pop_subcomp(EEG, final_remove, 0);
            EEG = eeg_checkset(EEG);
        end

        % 3. Average reference
        EEG = pop_reref(EEG, []);
        EEG = eeg_checkset(EEG);

        % 4. Save
        EEG.setname = [sid '_step2'];
        EEG = pop_saveset(EEG, 'filename', [sid '_step2.set'], 'filepath', paths.LAB_DATA);

        % 5. Over-correction QC vs step1
        flag_over = false;
        beta_pct = NaN;
        if exist(step1_path, 'file')
            try
                EEG_pre = pop_loadset('filename', [sid '_step1.set'], 'filepath', paths.LAB_DATA);
                qc = qc_ica_overcorrection(EEG_pre, EEG, sid);
                beta_pct = qc.beta_pct_change;
                flag_over = qc.flag;
                if flag_over
                    overcorr_flags{end+1} = sid; %#ok<AGROW>
                end
                % QC png
                try
                    qc_psd_plot(EEG_pre, EEG, sid, 'step2_ica_compare', paths.QC_PNG);
                catch
                end
            catch ME_qc
                fprintf('(over-corr QC fail: %s) ', ME_qc.message);
            end
        end

        % 6. Manifest
        info = struct();
        info.input            = ica_path;
        info.output           = out_path;
        info.final_remove     = final_remove;
        info.n_removed        = numel(final_remove);
        info.avg_reference    = true;
        info.beta_pct_change  = beta_pct;
        info.overcorr_flag    = flag_over;
        manifest_writer(sid, 'step2b', info);

        % 7. Log
        log_writer('step2b', sid, ...
            {'n_removed','remove_ic','beta_pct_change','overcorr_flag','status','note'}, ...
            {numel(final_remove), num2str(final_remove), beta_pct, flag_over, 'ok', ''});

        if flag_over
            fprintf('OK (removed %d, Beta %.1f%% ⚠️ OVER-CORR)\n', numel(final_remove), beta_pct);
        else
            fprintf('OK (removed %d, Beta %.1f%%)\n', numel(final_remove), beta_pct);
        end
        n_ok = n_ok + 1;
    catch ME
        log_writer('step2b', sid, ...
            {'n_removed','remove_ic','beta_pct_change','overcorr_flag','status','note'}, ...
            {NaN, '', NaN, false, 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

elapsed = toc(t_start);
fprintf('\n============================================================\n');
fprintf('  Step 2b done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Elapsed: %.1f min\n', elapsed/60);
fprintf('  Over-correction flags: %d  (%.1f%%) — threshold ±%d%%\n', ...
        numel(overcorr_flags), ...
        100*numel(overcorr_flags)/max(1,n_ok), thr.overcorr.warn_pct_change);
if ~isempty(overcorr_flags)
    fprintf('  Subjects: %s\n', strjoin(overcorr_flags, ', '));
    % 寫 warn list
    warn_path = fullfile(paths.LOGS, 'step2b_overcorr_warn.txt');
    fid = fopen(warn_path, 'w');
    fprintf(fid, '# Over-correction flags (Beta %% change > %d%%)\n', thr.overcorr.warn_pct_change);
    for k = 1:numel(overcorr_flags)
        fprintf(fid, '%s\n', overcorr_flags{k});
    end
    fclose(fid);
    fprintf('  Written: %s\n', warn_path);
end
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('============================================================\n');
end
