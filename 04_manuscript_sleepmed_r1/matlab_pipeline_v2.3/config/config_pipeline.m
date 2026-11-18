function pipe = config_pipeline()
% CONFIG_PIPELINE  pipeline 版本資訊與 ICA 演算法切換(v2.2)
%
% USAGE
%   p = config_pipeline();
%   p.version  % 'v2.2'
%
% v2.2 ─ 2026/5/27

pipe = struct();
pipe.version       = 'v2.3';
pipe.date_designed = '2026-11-16';   % v2.2 設計於 2026-05-27；v2.3 的更動見 CHANGELOG_v2.3.md

% 19-channel resting EEG montage(post-rename T3/4/5/6 → T7/8/P7/P8)
pipe.channels_19 = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz',...
                    'C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

% Pipeline orchestration
pipe.steps = { ...
    'step1_preprocess', ...
    'step2_ica', ...
    'review_ica_smart', ...
    'step2b_apply_ica', ...
    'step3_epoch_reject', ...
    'review_epochs_smart', ...
    'step4_export_psd', ...
    'master_qc_report', ...
    'master_methods_export' };

% 預先排除的個案：編號屬於受試者資料，不放在公開版。
% 本機請在 config/excluded_subjects.txt 一行寫一個編號（該檔不進版本控制）。
pipe.excluded_subjects = {};
f_excl = fullfile(fileparts(mfilename('fullpath')), 'excluded_subjects.txt');
if exist(f_excl, 'file')
    ids = strtrim(string(splitlines(fileread(f_excl))));
    pipe.excluded_subjects = cellstr(ids(strlength(ids) > 0))';
end

% Parallel computing
pipe.use_parfor   = true;
pipe.parfor_workers = 4;   % 視 CPU 調整

% Reproducibility 基準
pipe.global_seed = 42;

% Tool versions(實作時 autoexport 會 fill in)
pipe.tools = struct();
pipe.tools.eeglab  = '';   % 由 version_stamp 自動偵測
pipe.tools.iclabel = '';
pipe.tools.cleanline = '';
pipe.tools.firfilt = '';
pipe.tools.clean_rawdata = '';
end
