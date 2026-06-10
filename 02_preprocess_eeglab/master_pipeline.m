function master_pipeline(sid_list, opts)
% MASTER_PIPELINE  Fully automated EEG pipeline orchestrator.
%
% USAGE
%   master_pipeline();                          % 全部
%   master_pipeline({'sub01','sub02'});
%   master_pipeline([], struct('overwrite', true));
%
% opts:
%   .overwrite      = false  → 各 step 是否覆蓋既有檔
%   .skip_review    = true   → 跳過互動式 review(全自動)
%   .start_from     = 'step1' → 從哪一步開始
%   .stop_after     = 'step5' → 到哪一步停
%
% This copy is intended for new collection batches: it defaults to full
% automation while keeping the original preprocessing logic.

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end
if ~isfield(opts, 'skip_review'); opts.skip_review = true; end
if ~isfield(opts, 'start_from'); opts.start_from = 'step1'; end
if ~isfield(opts, 'stop_after'); opts.stop_after = 'step5'; end

paths = config_paths();
pipe  = config_pipeline();

% 確保 EEGLAB 啟動
if ~exist('pop_biosig', 'file')
    addpath(paths.EEGLAB_DIR);
    eeglab nogui;
end
addpath(genpath(paths.PIPELINE_ROOT));

steps = {'step1', 'step2', 'review_ica', 'step2b', 'step3', 'review_epochs', 'step4', 'step5'};
start_idx = find(strcmp(steps, opts.start_from));
stop_idx  = find(strcmp(steps, opts.stop_after));
if isempty(start_idx); start_idx = 1; end
if isempty(stop_idx); stop_idx = numel(steps); end

fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║   Generic Automated EEG Pipeline                              ║\n');
fprintf('║   Pipeline %s | %d subjects                                  ║\n', ...
        pipe.version, numel(sid_list));
fprintf('║   From: %-8s   To: %-8s   Skip review: %d              ║\n', ...
        opts.start_from, opts.stop_after, opts.skip_review);
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

stamp = version_stamp();
fprintf('Version stamp: %s\n', stamp.label);

t0 = tic;
for i = start_idx:stop_idx
    step = steps{i};
    fprintf('\n>>> Running %s ...\n', step);
    try
        switch step
            case 'step1'
                step1_preprocess(sid_list, opts);
            case 'step2'
                step2_ica(sid_list, opts);
            case 'review_ica'
                if opts.skip_review
                    fprintf('  (Skipped — skip_review=true; auto-confirm all ICLabel)\n');
                    auto_confirm_ica_decisions(sid_list, opts);
                else
                    review_ica_smart(sid_list);
                end
            case 'step2b'
                step2b_apply_ica(sid_list, opts);
            case 'step3'
                step3_epoch_reject(sid_list, opts);
            case 'review_epochs'
                if opts.skip_review
                    fprintf('  (Skipped — skip_review=true; epoch sampling already done in step3)\n');
                else
                    review_epochs_smart();
                end
            case 'step4'
                step4_export_psd(sid_list, opts);
            case 'step5'
                step5_export_bandpower(sid_list, opts);
        end
    catch ME
        fprintf('\n!!! Step %s failed: %s\n', step, ME.message);
        fprintf('Pipeline halted.\n');
        return;
    end
end

elapsed = toc(t0);
fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║   Master Pipeline DONE                                        ║\n');
fprintf('║   Total elapsed: %.1f min                                      ║\n', elapsed/60);
fprintf('║   Output root: %s\n', paths.PIPELINE_ROOT);
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
end


function auto_confirm_ica_decisions(sid_list, opts)
% Helper: 對所有 sid 自動寫 ica_decisions.json,直接套 ICLabel 自動 flag
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
for i = 1:numel(sid_list)
    sid = sid_list{i};
    ica_path = fullfile(paths.LAB_DATA, [sid '_ica.set']);
    dec_path = fullfile(paths.DECISIONS, [sid '_ica_decisions.json']);
    if ~exist(ica_path, 'file') || (exist(dec_path, 'file') && ~opts.overwrite)
        continue;
    end
    EEG = pop_loadset('filename', [sid '_ica.set'], 'filepath', paths.LAB_DATA);
    auto_flag = EEG.etc.iclabel_auto_flag;
    dec = struct();
    dec.subject        = sid;
    dec.timestamp      = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    dec.final_remove   = auto_flag;
    dec.auto_remove    = auto_flag;
    dec.auto_keep      = [];
    dec.review_queue   = [];
    dec.reviewer_notes = {'auto_confirm_no_review'};
    dec.version_stamp  = version_stamp();
    fid = fopen(dec_path, 'w');
    fprintf(fid, '%s', jsonencode(dec, 'PrettyPrint', true));
    fclose(fid);
end
end
