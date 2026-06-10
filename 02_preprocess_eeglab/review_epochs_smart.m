function review_epochs_smart(sid_list)
% REVIEW_EPOCHS_SMART  Smart epoch review(只對 low-quality 個案開 GUI)
%
% USAGE
%   review_epochs_smart();                  % 從 step3_need_review.txt 讀
%   review_epochs_smart({'sub03'});       % 特定個案
%
% 流程:
%   1. 讀 step3_need_review.txt(或傳入 list)
%   2. 對每人開 eegplot 視窗 → 點擊標紅 → 關閉確認
%   3. 若手動移除後仍 ≥ 60,重新採樣 60(per-subject seed,deterministic)
%   4. 覆蓋 step3.set + 更新 epoch_quality.json + log
%
% v2.2 ─ 2026/5/27

paths = config_paths();
thr   = config_thresholds();

% 取得 review queue
if nargin < 1 || isempty(sid_list)
    review_path = fullfile(paths.LOGS, 'step3_need_review.txt');
    if ~exist(review_path, 'file')
        fprintf('No review queue found (%s).\n', review_path);
        fprintf('Run step3_epoch_reject() first.\n');
        return;
    end
    sid_list = {};
    fid = fopen(review_path, 'r');
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line, '#'); continue; end
        % "sub03(20 epochs)" → 抽 sub03
        tok = regexp(line, '^(S\d+)', 'tokens', 'once');
        if ~isempty(tok)
            sid_list{end+1} = tok{1}; %#ok<AGROW>
        end
    end
    fclose(fid);
    sid_list = unique(sid_list);
end

fprintf('============================================================\n');
fprintf('  Smart Epoch Review\n');
fprintf('  Queue: %d subjects\n', numel(sid_list));
fprintf('  Click bad epochs (red) → close window → REJECT or CANCEL\n');
fprintf('============================================================\n\n');

for i = 1:numel(sid_list)
    sid = sid_list{i};
    set_path = fullfile(paths.LAB_DATA, [sid '_step3.set']);

    if ~exist(set_path, 'file')
        fprintf('[%d/%d] %s: no step3.set, skip\n', i, numel(sid_list), sid);
        continue;
    end

    fprintf('\n[%d/%d] %s\n', i, numel(sid_list), sid);

    try
        EEG = pop_loadset('filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);
        n_before = EEG.trials;
        fprintf('  Loaded %d epochs\n', n_before);

        % 開 eegplot
        DISPLAY_SCALE = 50;
        EPOCHS_PER_PAGE = 5;
        WIN_LENGTH = EPOCHS_PER_PAGE * thr.epoch.length_sec;
        eegplot(EEG.data, ...
            'srate',      EEG.srate, ...
            'title',      sprintf('%s — %d epochs', sid, n_before), ...
            'limits',     [EEG.times(1) EEG.times(end)], ...
            'command',    'EEG.reject.rejmanual = TMPREJ(:,1) > 0;', ...
            'eloc_file',  EEG.chanlocs, ...
            'winlength',  WIN_LENGTH, ...
            'dispchans',  19, ...
            'spacing',    DISPLAY_SCALE, ...
            'butlabel',   'REJECT marked');
        waitfor(gcf);

        n_removed = 0;
        if isfield(EEG.reject, 'rejmanual') && any(EEG.reject.rejmanual)
            rej_idx = find(EEG.reject.rejmanual);
            n_removed = length(rej_idx);
            EEG = pop_rejepoch(EEG, rej_idx, 0);
            EEG = eeg_checkset(EEG);
        end
        n_after = EEG.trials;

        % 重採樣 60
        seed = set_subject_rng(sid);
        target = thr.epoch.target_count;
        if n_after == 0
            error('All epochs removed during review');
        elseif n_after < target
            sel_idx = 1:n_after;
            fprintf('  Only %d remain (< %d), take all\n', n_after, target);
        else
            sel_idx = sort(randperm(n_after, target));
            fprintf('  Sampled %d from %d\n', target, n_after);
        end
        EEG = pop_select(EEG, 'trial', sel_idx);
        n_final = EEG.trials;

        EEG = pop_saveset(EEG, 'filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);

        % 更新 epoch_quality.json
        eq_path = fullfile(paths.DECISIONS, [sid '_epoch_quality.json']);
        if exist(eq_path, 'file')
            eq = jsondecode(fileread(eq_path));
        else
            eq = struct('subject', sid);
        end
        eq.review_n_before  = n_before;
        eq.review_n_removed = n_removed;
        eq.review_n_after   = n_after;
        eq.review_n_final   = n_final;
        eq.review_timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        eq.reviewed         = true;
        fid = fopen(eq_path, 'w');
        fprintf(fid, '%s', jsonencode(eq, 'PrettyPrint', true));
        fclose(fid);

        log_writer('review_epochs', sid, ...
            {'n_before','n_removed','n_after','n_final','status','note'}, ...
            {n_before, n_removed, n_after, n_final, 'ok', ''});

        fprintf('  Reviewed: %d → removed %d → %d remain → final %d\n', ...
                n_before, n_removed, n_after, n_final);
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        log_writer('review_epochs', sid, ...
            {'n_before','n_removed','n_after','n_final','status','note'}, ...
            {NaN, NaN, NaN, NaN, 'error', ME.message});
    end
end

close all;
fprintf('\n=== Epoch review done. Next: step4_export_psd() ===\n');
end
