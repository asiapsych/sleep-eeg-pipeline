function pipe = config_pipeline()
% CONFIG_PIPELINE  Pipeline metadata and behavior switches.
%
% USAGE
%   p = config_pipeline();
%   p.version  % 'v2.2'
%
pipe = struct();
pipe.version       = 'generic-auto-v1.0';
pipe.date_designed = '2026-06-03';

% 19-channel resting EEG montage(post-rename T3/4/5/6 → T7/8/P7/P8)
pipe.channels_19 = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz',...
                    'C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

% Pipeline orchestration.
pipe.steps = { ...
    'step1_preprocess', ...
    'step2_ica', ...
    'step2b_apply_ica', ...
    'step3_epoch_reject', ...
    'step4_export_psd', ...
    'step5_export_bandpower', ...
    'master_qc_report', ...
    'master_methods_export' };

% Per-project exclusions. Keep empty for new collection unless a subject is invalid.
pipe.excluded_subjects = {};

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
