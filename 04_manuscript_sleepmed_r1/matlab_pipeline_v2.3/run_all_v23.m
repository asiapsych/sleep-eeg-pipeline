% run_all_v23.m  v2.3 完整流程：原始 EDF → step1 → step2 → step2b → step3 → step4 → step5（全自動、固定門檻）
%
% 執行：
%   matlab -nodisplay -nosplash -nodesktop -batch "run('run_all_v23.m')"
% 輸出：本資料夾 outputs/（Lab_data、decisions、manifests、logs、psd_mat）
% 紀錄：outputs/logs/run_all_v23_<時間>.log（完整輸出）、各步驟的 step*_log.csv、每人的 manifests/<SID>.json
% 不畫圖、不呼叫互動式審查。任何一步對某位受試者失敗，都會寫進紀錄，不會靜靜改用別的方法。
% 2026-11-16

here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
p0 = config_paths(); addpath(p0.EEGLAB_DIR);
eeglab nogui;
p = config_paths(); pipe = config_pipeline(); thr = config_thresholds();
assert(contains(p.PIPELINE_ROOT, 'matlab_pipeline_v2.3'), 'wrong pipeline root');
assert(contains(which('step2_ica'), 'matlab_pipeline_v2.3'), 'wrong step2_ica on path');
assert(strcmp(pipe.version, 'v2.3') && strcmpi(thr.ica.algorithm, 'picard'));

if exist('V23_TEST_SIDS', 'var'); sids = V23_TEST_SIDS; else; sids = subjects(); end
stamp = version_stamp();
fprintf('v2.3 full run: %d subjects, %s\n%s\n', numel(sids), datestr(now), stamp.label);
disp(stamp); disp(stamp.plugins);
t0 = tic;
master_pipeline(sids, struct('overwrite', false, 'skip_review', true, 'start_from', 'step1', 'stop_after', 'step4'));
step5_export_epoch_psd(sids);   % 逐 epoch 頻譜（清醒度的逐段分析用）
fprintf('v2.3 full run done: %.1f min, %s\n', toc(t0) / 60, datestr(now));
