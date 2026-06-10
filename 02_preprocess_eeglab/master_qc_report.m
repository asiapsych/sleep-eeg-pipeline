function master_qc_report(sid_list, opts)
% MASTER_QC_REPORT  Per-subject multi-page QC PDF report
%
% USAGE
%   master_qc_report();
%   master_qc_report({'sub01'});
%
% 對每個 subject 產一份 PDF(qc_pdf/<sid>.pdf),包含:
%   Page 1:Header + summary(version stamp + grade + 主要參數)
%   Page 2:Step 1 PSD before/after + bad channel topography
%   Page 3:Step 2 ICA 各 IC 之 ICLabel scores + auto-decision
%   Page 4:Step 2b over-correction QC(per-band % change)
%   Page 5:Step 3 epoch rejection waterfall + sampled epoch PSD
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
thr   = config_thresholds();
bands = config_bands();

fprintf('============================================================\n');
fprintf('  Master QC Report\n');
fprintf('  Generating PDF per subject (%d total)\n', numel(sid_list));
fprintf('============================================================\n\n');

n_ok = 0; n_skip = 0; n_fail = 0;
for i = 1:numel(sid_list)
    sid = sid_list{i};
    pdf_path = fullfile(paths.QC_PDF, [sid '_qc.pdf']);

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if exist(pdf_path, 'file') && ~opts.overwrite
        fprintf('SKIP (done)\n');
        n_skip = n_skip + 1;
        continue;
    end

    mf_path = fullfile(paths.MANIFEST_DIR, [sid '.json']);
    if ~exist(mf_path, 'file')
        fprintf('SKIP (no manifest)\n');
        n_skip = n_skip + 1;
        continue;
    end

    try
        mf = jsondecode(fileread(mf_path));

        % --- Page 1: Header / Summary ---
        fig = figure('Visible', 'off', 'Position', [0 0 850 1100], 'Color', 'w');
        axis off;
        grade = compute_grade(mf, thr);
        txt_lines = {
            sprintf('SUBJECT QC REPORT — %s', sid)
            sprintf('Pipeline %s', getfield_safe(mf, {'version_stamp', 'label'}, '?'))
            ''
            sprintf('Quality Grade: %s', grade)
            ''
            '--- STEP 1 (Load + Bandpass + Bad Channel) ---'
            sprintf('  n_chan_in: %d, srate: %g Hz, duration: %.1f s', ...
                getfield_safe(mf, {'step1','n_chan_in'}, NaN), ...
                getfield_safe(mf, {'step1','srate'}, NaN), ...
                getfield_safe(mf, {'step1','duration_s'}, NaN))
            sprintf('  bandpass: %s', mat2str(getfield_safe(mf, {'step1','bandpass'}, [NaN NaN])))
            sprintf('  bad channels (n=%d): %s', ...
                getfield_safe(mf, {'step1','n_bad'}, 0), ...
                safe_join_strings(getfield_safe(mf, {'step1','bad_channels'}, {})))
            ''
            '--- STEP 2 (CleanLine + ICA + ICLabel) ---'
            sprintf('  ICA algorithm: %s, rank: %d, components: %d', ...
                getfield_safe(mf, {'step2','ica_algorithm'}, '?'), ...
                getfield_safe(mf, {'step2','data_rank'}, NaN), ...
                getfield_safe(mf, {'step2','n_components'}, NaN))
            sprintf('  ICLabel auto-flag: %d', getfield_safe(mf, {'step2','iclabel_flag_count'}, 0))
            ''
            '--- STEP 2b (Apply + Avg Ref + Over-correction QC) ---'
            sprintf('  IC removed: %d', getfield_safe(mf, {'step2b','n_removed'}, 0))
            sprintf('  Beta %% change (before/after ICA): %.1f%%', ...
                getfield_safe(mf, {'step2b','beta_pct_change'}, NaN))
            sprintf('  Over-correction flag: %d (threshold ±%d%%)', ...
                getfield_safe(mf, {'step2b','overcorr_flag'}, 0), thr.overcorr.warn_pct_change)
            ''
            '--- STEP 3 (Epoch + 3-layer reject + sample 60) ---'
            sprintf('  epochs: total=%d → A=%d → B=%d → C=%d → final=%d', ...
                getfield_safe(mf, {'step3','n_total'}, NaN), ...
                getfield_safe(mf, {'step3','n_after_layerA'}, NaN), ...
                getfield_safe(mf, {'step3','n_after_layerB'}, NaN), ...
                getfield_safe(mf, {'step3','n_after_layerC'}, NaN), ...
                getfield_safe(mf, {'step3','n_final'}, NaN))
            sprintf('  reject rate: %.0f%%, below target: %d', ...
                getfield_safe(mf, {'step3','reject_rate'}, NaN)*100, ...
                getfield_safe(mf, {'step3','below_target'}, 0))
            ''
            '--- STEP 4 (PSD export) ---'
            sprintf('  PSD shape: [n_ch=%d × n_freq=%d], n_ep=%d', ...
                19, ...
                getfield_safe(mf, {'step4','n_freqs'}, NaN), ...
                getfield_safe(mf, {'step4','n_ep'}, NaN))
            ''
            sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'))
        };
        text(0.05, 0.95, txt_lines, 'FontName', 'Courier', 'FontSize', 9, ...
             'VerticalAlignment', 'top', 'Interpreter', 'none');
        title(sprintf('QC Report — %s', sid), 'FontSize', 16);
        exportgraphics(fig, pdf_path, 'ContentType', 'vector', 'Append', false);

        % --- Page 2: Step 1 PSD before/after + bad ch topography ---
        step1_png = fullfile(paths.QC_PNG, sprintf('%s_step1_psd_compare.png', sid));
        if exist(step1_png, 'file')
            clf(fig);
            img = imread(step1_png);
            imshow(img, 'Border', 'tight');
            title(sprintf('%s — Step 1 (Bandpass)', sid));
            exportgraphics(fig, pdf_path, 'Append', true);
        end

        % --- Page 3: Step 2b ICA over-correction ---
        step2_png = fullfile(paths.QC_PNG, sprintf('%s_step2_ica_compare_psd_compare.png', sid));
        if exist(step2_png, 'file')
            clf(fig);
            img = imread(step2_png);
            imshow(img, 'Border', 'tight');
            title(sprintf('%s — Step 2 (ICA before/after)', sid));
            exportgraphics(fig, pdf_path, 'Append', true);
        end

        % --- Page 4: Epoch rejection waterfall ---
        if isfield(mf, 'step3')
            s3 = mf.step3;
            clf(fig);
            labels_e = {'Total', 'After A', 'After B', 'After C', 'Final'};
            counts_e = [getfield_safe(s3,{'n_total'},0), ...
                        getfield_safe(s3,{'n_after_layerA'},0), ...
                        getfield_safe(s3,{'n_after_layerB'},0), ...
                        getfield_safe(s3,{'n_after_layerC'},0), ...
                        getfield_safe(s3,{'n_final'},0)];
            bar(counts_e);
            set(gca, 'XTickLabel', labels_e);
            ylabel('# epochs');
            title(sprintf('%s — Epoch Rejection Waterfall', sid));
            grid on;
            exportgraphics(fig, pdf_path, 'Append', true);
        end

        % --- Page 5: Final step3 PSD ---
        step3_set = fullfile(paths.LAB_DATA, [sid '_step3.set']);
        if exist(step3_set, 'file')
            try
                EEG = pop_loadset('filename', [sid '_step3.set'], 'filepath', paths.LAB_DATA);
                clf(fig);
                data = reshape(EEG.data, EEG.nbchan, []);
                psd_all = [];
                for ch = 1:EEG.nbchan
                    [pxx, ff] = pwelch(data(ch,:), hamming(thr.psd.win_samples), 0, ...
                                        thr.psd.nfft, EEG.srate);
                    if isempty(psd_all); psd_all = zeros(EEG.nbchan, length(pxx)); end
                    psd_all(ch,:) = pxx;
                end
                m = (ff >= 1) & (ff <= 50);
                semilogy(ff(m), psd_all(:,m)', 'LineWidth', 0.7);
                hold on;
                semilogy(ff(m), mean(psd_all(:,m),1), 'k-', 'LineWidth', 2.5);
                xlabel('Frequency (Hz)'); ylabel('PSD (\muV^2/Hz)');
                title(sprintf('%s — Final step3 PSD (per ch + mean)', sid));
                xline(4,'--','D'); xline(8,'--','T'); xline(13,'--','A');
                xline(30,'--','B'); xline(45,'--','G');
                grid on; legend([{EEG.chanlocs.labels}, {'Mean'}], 'Location', 'eastoutside', 'FontSize', 7);
                exportgraphics(fig, pdf_path, 'Append', true);
            catch
            end
        end

        close(fig);
        fprintf('OK\n');
        n_ok = n_ok + 1;
    catch ME
        fprintf('FAIL: %s\n', ME.message);
        if exist('fig', 'var') && ishandle(fig); close(fig); end
        n_fail = n_fail + 1;
    end
end

fprintf('\nQC Report done: %d ok, %d skip, %d fail / %d total\n', n_ok, n_skip, n_fail, numel(sid_list));
fprintf('Output: %s\n', paths.QC_PDF);
end


function grade = compute_grade(mf, thr)
% A/B/C/F grade
g = 'A';
n_bad = getfield_safe(mf, {'step1','n_bad'}, 0);
n_ic_removed = getfield_safe(mf, {'step2b','n_removed'}, 0);
n_final = getfield_safe(mf, {'step3','n_final'}, 0);
over_flag = getfield_safe(mf, {'step2b','overcorr_flag'}, 0);

if n_bad > thr.qc.grade_B_max_bad_ch || ...
   n_ic_removed > thr.qc.grade_B_max_ic_removed || ...
   n_final < thr.qc.grade_B_min_epoch
    g = 'F';
elseif n_bad > thr.qc.grade_A_max_bad_ch || ...
       n_ic_removed > thr.qc.grade_A_max_ic_removed || ...
       n_final < thr.qc.grade_A_min_epoch
    g = 'C';
elseif over_flag
    g = 'B';  % Over-correction = B
end
grade = g;
end


function v = getfield_safe(s, path, default)
% Safely traverse nested fields
v = default;
try
    for k = 1:numel(path)
        s = s.(path{k});
    end
    v = s;
catch
end
end


function s = safe_join_strings(x)
% 處理各種型別,安全轉成 ', ' 分隔字串
if isempty(x)
    s = '(none)';
elseif iscell(x)
    parts = {};
    for k = 1:numel(x)
        if ischar(x{k}) || isstring(x{k})
            parts{end+1} = char(x{k}); %#ok<AGROW>
        elseif isnumeric(x{k})
            parts{end+1} = num2str(x{k}); %#ok<AGROW>
        end
    end
    if isempty(parts); s = '(none)'; else; s = strjoin(parts, ', '); end
elseif ischar(x)
    s = x;
elseif isstring(x)
    s = char(strjoin(x, ', '));
elseif isnumeric(x)
    s = num2str(x(:)');
else
    s = '(unknown)';
end
end
