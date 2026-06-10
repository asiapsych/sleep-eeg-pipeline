function fig = qc_psd_plot(EEG_before, EEG_after, sid, label, out_path)
% QC_PSD_PLOT  畫 PSD before/after 比對(per channel + per band % change)
%
% USAGE
%   qc_psd_plot(EEG_pre, EEG_post, 'sub01', 'step2_ica', '/path/to/qc_png/');
%
% 輸出:
%   {out_path}/{sid}_{label}_psd_compare.png
%
% v2.2 ─ 2026/5/27

bands = config_bands();

% 計算兩個 EEG 之 PSD
[psd_b, f_b] = compute_psd_brain_avg(EEG_before);
[psd_a, f_a] = compute_psd_brain_avg(EEG_after);

% 對齊 freq 範圍
f_min = max(min(f_b), min(f_a));
f_max = min(max(f_b), max(f_a));
m_b = (f_b >= f_min) & (f_b <= f_max);
m_a = (f_a >= f_min) & (f_a <= f_max);
psd_b = psd_b(m_b); f_b = f_b(m_b);
psd_a = psd_a(m_a); f_a = f_a(m_a);

% Figure
fig = figure('Visible', 'off', 'Position', [100 100 1200 700]);

% Panel 1: PSD 對比
subplot(2, 2, [1 2]);
semilogy(f_b, psd_b, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Before'); hold on;
semilogy(f_a, psd_a, 'r-', 'LineWidth', 1.5, 'DisplayName', 'After');
xline(4, '--', 'D');
xline(8, '--', 'T');
xline(13, '--', 'A');
xline(30, '--', 'B');
xline(45, '--', 'G');
xlim([1 70]);
xlabel('Frequency (Hz)');
ylabel('PSD (\muV^2/Hz)');
title(sprintf('%s — %s PSD (BRAIN avg)', sid, label));
legend('Location', 'best');
grid on;

% Panel 2: % change per primary band
subplot(2, 2, 3);
primary_names = fieldnames(bands.primary);
n_band = numel(primary_names);
pct_change = zeros(1, n_band);
band_powers_b = zeros(1, n_band);
band_powers_a = zeros(1, n_band);
for i = 1:n_band
    rng = bands.primary.(primary_names{i});
    pb = trapz(f_b((f_b>=rng(1))&(f_b<=rng(2))), psd_b((f_b>=rng(1))&(f_b<=rng(2))));
    pa = trapz(f_a((f_a>=rng(1))&(f_a<=rng(2))), psd_a((f_a>=rng(1))&(f_a<=rng(2))));
    band_powers_b(i) = pb;
    band_powers_a(i) = pa;
    if pb > 0
        pct_change(i) = (pa - pb) / pb * 100;
    end
end
bar_colors = repmat([0.4 0.4 0.7], n_band, 1);
% Highlight Beta band
beta_idx = find(strcmp(primary_names, 'Beta'));
if ~isempty(beta_idx) && abs(pct_change(beta_idx)) > 20
    bar_colors(beta_idx, :) = [0.8 0.2 0.2];  % 紅旗
end
b = bar(pct_change, 'FaceColor', 'flat');
b.CData = bar_colors;
set(gca, 'XTickLabel', primary_names);
ylabel('% change (After - Before) / Before');
title('Per-band % change');
grid on;
yline(0, '-k');
yline(20, '--r');
yline(-20, '--r');

% Panel 3: Abs power(linear)
subplot(2, 2, 4);
bar([band_powers_b; band_powers_a]');
set(gca, 'XTickLabel', primary_names);
ylabel('Absolute power (\muV^2)');
legend({'Before', 'After'}, 'Location', 'best');
title('Absolute band power');
grid on;

% 存檔
if ~isempty(out_path)
    if ~exist(out_path, 'dir'); mkdir(out_path); end
    out_file = fullfile(out_path, sprintf('%s_%s_psd_compare.png', sid, label));
    saveas(fig, out_file);
end
end


function [psd, f] = compute_psd_brain_avg(EEG)
% 對連續訊號做 pwelch,跨 channel 平均
thr = config_thresholds();
data = EEG.data;
if ndims(data) == 3  % epoched
    data = reshape(data, EEG.nbchan, []);
end
nperseg = min(thr.psd.win_samples, size(data, 2));
psd_all = zeros(EEG.nbchan, floor(nperseg/2)+1);
for ci = 1:EEG.nbchan
    [pxx, f] = pwelch(data(ci, :), hamming(nperseg), 0, nperseg, EEG.srate);
    if ci == 1
        psd_all = zeros(EEG.nbchan, length(pxx));
    end
    psd_all(ci, :) = pxx';
end
psd = mean(psd_all, 1);
end
