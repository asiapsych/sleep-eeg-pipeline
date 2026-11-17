function validation_runica_vs_picard(sample_sids, opts)
% VALIDATION_RUNICA_VS_PICARD  Sensitivity:對 N=20 個 subject 跑 runica + picard 各 1 次,比較
%
% USAGE
%   validation_runica_vs_picard();                  % 隨機抽 20 人
%   validation_runica_vs_picard({'S000001', ...}, ...);
%
% 對每 subject:
%   1. 跑 runica + ICLabel
%   2. 跑 picard + ICLabel
%   3. 比較:IC topography(每 IC 對最相似的 r);ICLabel scores 相關性
%
% 輸出:
%   outputs/methods/validation_runica_vs_picard.csv
%   outputs/methods/validation_runica_vs_picard.png
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sample_sids)
    rng(42, 'twister');
    all_sids = subjects();
    sample_sids = all_sids(randperm(numel(all_sids), min(20, numel(all_sids))));
end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'check_picard'); opts.check_picard = true; end

paths = config_paths();
thr   = config_thresholds();

% 檢查 picard 是否裝
has_picard = ~isempty(which('picard'));
if ~has_picard
    fprintf('\n[警告] picard plugin 未裝(在 EEGLAB GUI → File → Manage Extensions 安裝)\n');
    fprintf('此 validation 需 picard。先 skip。\n');
    return;
end

fprintf('============================================================\n');
fprintf('  Sensitivity:runica vs picard\n');
fprintf('  Sample subjects: %d\n', numel(sample_sids));
fprintf('============================================================\n\n');

results = struct('subject', {}, 'topo_r_mean', {}, 'topo_r_min', {}, ...
                 'iclabel_brain_r', {});

for i = 1:numel(sample_sids)
    sid = sample_sids{i};
    in_path = fullfile(paths.LAB_DATA, [sid '_step1.set']);
    if ~exist(in_path, 'file')
        fprintf('[%d/%d] %s: no step1.set, skip\n', i, numel(sample_sids), sid);
        continue;
    end

    fprintf('[%d/%d] %s ... ', i, numel(sample_sids), sid);

    try
        % Load + CleanLine
        EEG = pop_loadset('filename', [sid '_step1.set'], 'filepath', paths.LAB_DATA);
        EEG = pop_cleanline(EEG, ...
            'bandwidth', thr.cleanline.bandwidth, ...
            'chanlist', 1:EEG.nbchan, ...
            'linefreqs', thr.cleanline.line_freq, ...
            'p', 0.01, 'pad', 2, 'plotfigures', 0, ...
            'scanforlines', thr.cleanline.scanForLines, ...
            'sigtype', 'Channels', 'tau', thr.cleanline.tau, ...
            'verb', 0, 'winsize', 4, 'winstep', 1, ...
            'normSpectrum', 0, 'computepower', 0);

        % 估 rank
        mf_path = fullfile(paths.MANIFEST_DIR, [sid '.json']);
        n_interp = 0;
        if exist(mf_path, 'file')
            try
                mf = jsondecode(fileread(mf_path));
                n_interp = safef(mf, {'step1','n_bad'}, 0);
            catch
            end
        end
        data_rank = compute_data_rank(EEG, n_interp);

        % runica
        seed = set_subject_rng(sid);
        EEG_r = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, ...
                           'pca', data_rank, 'interupt', 'off');
        EEG_r = iclabel(EEG_r);

        % picard(EEGLAB 標準參數)
        set_subject_rng(sid);
        EEG_p = pop_runica(EEG, 'icatype', 'picard', 'maxiter', 500, ...
                           'mode', 'standard', 'pca', data_rank, 'interupt', 'off');
        EEG_p = iclabel(EEG_p);

        % 比 topography:對每 runica IC,找 picard IC 最相似的(max |r|)
        winv_r = EEG_r.icawinv;   % nbchan × n_comp
        winv_p = EEG_p.icawinv;
        n_comp = size(winv_r, 2);
        topo_rs = zeros(1, n_comp);
        for c = 1:n_comp
            rs = abs(corr(winv_r(:,c), winv_p));
            topo_rs(c) = max(rs);
        end

        % 比 ICLabel Brain scores 之相關
        cls_r = EEG_r.etc.ic_classification.ICLabel.classes;
        sc_r = EEG_r.etc.ic_classification.ICLabel.classifications;
        sc_p = EEG_p.etc.ic_classification.ICLabel.classifications;
        brain_col = find(strcmpi(cls_r, 'Brain'));
        % 對應 IC 順序(不一定一樣)— 用 topo best match
        brain_r = NaN;
        if ~isempty(brain_col)
            sc_p_matched = zeros(n_comp, 1);
            for c = 1:n_comp
                rs = abs(corr(winv_r(:,c), winv_p));
                [~, best] = max(rs);
                sc_p_matched(c) = sc_p(best, brain_col);
            end
            brain_r = corr(sc_r(:, brain_col), sc_p_matched);
        end

        results(end+1) = struct('subject', sid, ...
            'topo_r_mean', mean(topo_rs), ...
            'topo_r_min', min(topo_rs), ...
            'iclabel_brain_r', brain_r); %#ok<AGROW>

        fprintf('topo_r mean=%.3f min=%.3f Brain_r=%.3f\n', ...
                mean(topo_rs), min(topo_rs), brain_r);
    catch ME
        fprintf('FAIL: %s\n', ME.message);
    end
end

% 寫 CSV
out_csv = fullfile(paths.METHODS_OUT, 'validation_runica_vs_picard.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, 'subject,topo_r_mean,topo_r_min,iclabel_brain_r\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '%s,%.4f,%.4f,%.4f\n', r.subject, r.topo_r_mean, r.topo_r_min, r.iclabel_brain_r);
end
fclose(fid);
fprintf('\nResults: %s\n', out_csv);

% Summary
if ~isempty(results)
    fprintf('\n--- Summary (N=%d) ---\n', numel(results));
    fprintf('Topo r (mean across IC): mean=%.3f, SD=%.3f\n', ...
            mean([results.topo_r_mean]), std([results.topo_r_mean]));
    fprintf('ICLabel Brain r (subject-level): mean=%.3f, SD=%.3f\n', ...
            nanmean([results.iclabel_brain_r]), nanstd([results.iclabel_brain_r]));
end
end


function v = safef(s, path, default)
v = default;
try
    for k = 1:numel(path); s = s.(path{k}); end
    v = double(s);
catch
end
end
