function validation_test_retest()
% VALIDATION_TEST_RETEST  T1/T2/T3 Beta power ICC reliability
%
% USAGE
%   validation_test_retest();
%
% 對「至少 2 個 timepoint 有 EEG」之個案,算 Beta power 之 ICC(2,1):
%   - 跨 channel 之 BRAIN-avg Beta_rel
%   - Test–retest reliability 在 T1-T2、T1-T3、T2-T3 各算
%
% 輸出:
%   outputs/methods/validation_test_retest.csv
%
% v2.2 ─ 2026/5/27
%
% 注意:目前版本只算 T1 EEG(因為我們現在只跑 T1 pipeline)。
% T2/T3 EEG 若有資料,需先跑 v2.2 pipeline(等 T2/T3 .ng 收集完成)。

paths = config_paths();
bands = config_bands();

sid_list = subjects();

% 從 PSD .mat 算每人 Beta_rel BRAIN avg
beta_rel = struct();
for i = 1:numel(sid_list)
    sid = sid_list{i};
    psd_path = fullfile(paths.PSD_OUT, sprintf('psd_%s.mat', sid));
    if ~exist(psd_path, 'file'); continue; end
    try
        d = load(psd_path);
        % per-channel total = sum of 5 PRIMARY bands
        primary_names = fieldnames(bands.primary);
        n_ch = size(d.psd, 1);
        ch_beta = zeros(1, n_ch);
        ch_total = zeros(1, n_ch);
        for ch = 1:n_ch
            for b = 1:numel(primary_names)
                rng = bands.primary.(primary_names{b});
                mask = (d.freqs >= rng(1)) & (d.freqs <= rng(2));
                bp = trapz(d.freqs(mask), d.psd(ch, mask));
                ch_total(ch) = ch_total(ch) + bp;
                if strcmp(primary_names{b}, 'Beta')
                    ch_beta(ch) = bp;
                end
            end
        end
        % BRAIN-avg relative power (pooled)
        beta_rel.(sid) = mean(ch_beta) / mean(ch_total) * 100;
    catch
    end
end

fprintf('Subjects with Beta_rel computed: %d\n', numel(fieldnames(beta_rel)));
fprintf('\nT2/T3 EEG pipeline 尚未跑完整時,test-retest ICC 無法直接算。\n');
fprintf('預留接口:當 T2/T3 之 psd_<sid>_T2.mat / _T3.mat 出現時,改 join 三時點即可。\n');
fprintf('\n目前先輸出 T1 Beta_rel CSV 供日後合併:\n');

out_csv = fullfile(paths.METHODS_OUT, 'validation_T1_beta_rel.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, 'NEW_NUMBER,Beta_rel_BRAIN_T1\n');
names = fieldnames(beta_rel);
for i = 1:numel(names)
    fprintf(fid, '%s,%.4f\n', names{i}, beta_rel.(names{i}));
end
fclose(fid);
fprintf('  %s\n', out_csv);
end
