% run_step6_v23.m  只跑 step6（表面拉普拉斯後的頻譜；肌電的敏感度分析用）。前提：run_all_v23.m 已完成。
%
% 執行：
%   matlab -nodisplay -nosplash -nodesktop -batch "run('run_step6_v23.m')"
% 輸出：outputs/psd_mat_laplacian/psd_<SID>.mat；紀錄：outputs/logs/run_step6_v23_<時間>.log、step6_log.csv
% 不畫圖。2026-11-17

here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
p0 = config_paths(); addpath(p0.EEGLAB_DIR);
eeglab nogui;
p = config_paths();
assert(contains(p.PIPELINE_ROOT, 'matlab_pipeline_v2.3'), 'wrong pipeline root');
assert(contains(which('step6_export_psd_laplacian'), 'matlab_pipeline_v2.3'), 'wrong step6 on path');
diary(fullfile(p.LOGS, ['run_step6_v23_' datestr(now, 'yyyymmdd_HHMMSS') '.log']));
stamp = version_stamp(); fprintf('v2.3 step6: %s\n%s\n', datestr(now), stamp.label);
t0 = tic;
step6_export_psd_laplacian(subjects());
fprintf('v2.3 step6 done: %.1f min, %s\n', toc(t0) / 60, datestr(now));
diary off;
