function population_qc_summary()
% POPULATION_QC_SUMMARY  Single-PNG population overview
%
% USAGE
%   population_qc_summary();
%
% 讀所有 manifest,產 1 張 PNG 涵蓋:
%   - Pass rate per step
%   - Bad channel 分布直方圖
%   - ICA components removed 分布
%   - Epoch retention 分布
%   - Beta % change 分布
%   - Quality grade 分布
%
% v2.2 ─ 2026/5/27

paths = config_paths();
thr   = config_thresholds();

sid_list = subjects();
n_total = numel(sid_list);

agg = struct('n_bad', [], 'n_ic_removed', [], 'beta_pct', [], 'n_final', [], ...
             'overcorr', [], 'reject_rate', [], 'grade', '');
for i = 1:n_total
    mf_path = fullfile(paths.MANIFEST_DIR, [sid_list{i} '.json']);
    if ~exist(mf_path, 'file'); continue; end
    try
        mf = jsondecode(fileread(mf_path));
        agg.n_bad(end+1)        = safef(mf, {'step1','n_bad'}, NaN);
        agg.n_ic_removed(end+1) = safef(mf, {'step2b','n_removed'}, NaN);
        agg.beta_pct(end+1)     = safef(mf, {'step2b','beta_pct_change'}, NaN);
        agg.n_final(end+1)      = safef(mf, {'step3','n_final'}, NaN);
        agg.overcorr(end+1)     = safef(mf, {'step2b','overcorr_flag'}, 0);
        agg.reject_rate(end+1)  = safef(mf, {'step3','reject_rate'}, NaN);
    catch
    end
end
n_done = numel(agg.n_bad);

fig = figure('Visible', 'off', 'Position', [0 0 1400 900], 'Color', 'w');

% Panel 1: Pass rate per step
subplot(2, 3, 1);
n_step1 = sum(~isnan(agg.n_bad));
n_step2b = sum(~isnan(agg.n_ic_removed));
n_step3 = sum(~isnan(agg.n_final));
labels_st = {'Step1', 'Step2b', 'Step3'};
vals = [n_step1, n_step2b, n_step3];
bar(vals);
set(gca, 'XTickLabel', labels_st);
ylabel('# subjects completed');
title(sprintf('Pass rate per step (target N=%d)', n_total));
for k = 1:numel(vals)
    text(k, vals(k)+1, sprintf('%d', vals(k)), 'HorizontalAlignment', 'center');
end
ylim([0, n_total+5]);
grid on;

% Panel 2: Bad channel distribution
subplot(2, 3, 2);
histogram(agg.n_bad(~isnan(agg.n_bad)), 0:1:max(agg.n_bad)+1);
xlabel('# bad channels'); ylabel('# subjects');
title(sprintf('Bad channels (mean %.2f, SD %.2f)', nanmean(agg.n_bad), nanstd(agg.n_bad)));
grid on;

% Panel 3: ICA components removed
subplot(2, 3, 3);
histogram(agg.n_ic_removed(~isnan(agg.n_ic_removed)));
xlabel('# components removed'); ylabel('# subjects');
title(sprintf('ICA removed (mean %.2f, SD %.2f)', nanmean(agg.n_ic_removed), nanstd(agg.n_ic_removed)));
grid on;

% Panel 4: Epoch retention
subplot(2, 3, 4);
histogram(agg.n_final(~isnan(agg.n_final)));
xlabel('# final epochs'); ylabel('# subjects');
title(sprintf('Final epoch count (target %d)', thr.epoch.target_count));
xline(thr.epoch.target_count, '--r', 'target');
grid on;

% Panel 5: Beta % change distribution (over-correction)
subplot(2, 3, 5);
histogram(agg.beta_pct(~isnan(agg.beta_pct)));
xlabel('Beta % change (after-before ICA)'); ylabel('# subjects');
title(sprintf('Over-correction (red=flagged, %d of %d)', sum(agg.overcorr==1), n_done));
xline(thr.overcorr.warn_pct_change, '--r');
xline(-thr.overcorr.warn_pct_change, '--r');
grid on;

% Panel 6: Reject rate
subplot(2, 3, 6);
histogram(agg.reject_rate(~isnan(agg.reject_rate))*100);
xlabel('Reject rate %'); ylabel('# subjects');
title('Step3 epoch reject rate');
xline(thr.epoch.review_max_reject*100, '--r');
grid on;

sgtitle(sprintf('Population QC Summary — Pipeline v2.2 (N=%d/%d)', n_done, n_total), ...
        'FontSize', 14, 'FontWeight', 'bold');

out_path = fullfile(paths.QC_PNG, 'population_qc_summary.png');
saveas(fig, out_path);
close(fig);

fprintf('Population QC summary: %s\n', out_path);

% 也產 outliers top-10 worst
fprintf('\n--- Top 10 worst (most epochs rejected) ---\n');
[~, idx] = sort(agg.reject_rate, 'descend', 'MissingPlacement', 'last');
top_n = min(10, sum(~isnan(agg.reject_rate)));
sid_valid = sid_list(~isnan(agg.reject_rate));
agg_valid_rej = agg.reject_rate(~isnan(agg.reject_rate));
[sorted_rej, sorted_idx] = sort(agg_valid_rej, 'descend');
fprintf('%-12s  %-7s\n', 'subject', 'reject%');
for k = 1:top_n
    fprintf('%-12s  %.1f%%\n', sid_valid{sorted_idx(k)}, sorted_rej(k)*100);
end

% 也存 outliers CSV
csv_path = fullfile(paths.LOGS, 'population_outliers_top10.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'rank,subject,reject_pct\n');
for k = 1:top_n
    fprintf(fid, '%d,%s,%.3f\n', k, sid_valid{sorted_idx(k)}, sorted_rej(k));
end
fclose(fid);
fprintf('\nOutliers: %s\n', csv_path);
end


function v = safef(s, path, default)
v = default;
try
    for k = 1:numel(path)
        s = s.(path{k});
    end
    v = double(s);
catch
end
end
