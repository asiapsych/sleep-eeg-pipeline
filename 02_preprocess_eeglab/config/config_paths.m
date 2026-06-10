function paths = config_paths()
% CONFIG_PATHS  Centralized file paths for the generic automated EEG pipeline.
%
% USAGE
%   p = config_paths();
%   raw_dir = p.RAW_EDF_DIR;
%
% Edit config_user.m for each new collection/project instead of changing
% processing scripts.

paths = struct();

this_file = mfilename('fullpath');
paths.PIPELINE_ROOT = fileparts(fileparts(this_file));

% Input: put new EDF files here by default, or override in config_user.m.
paths.RAW_EDF_DIR = fullfile(paths.PIPELINE_ROOT, 'raw_edf');
paths.RAW_FILE_PATTERN = '*.edf';
paths.RAW_NG_DIR  = '';
paths.NG_MANIFEST = '';

% Working/output folders.
paths.LAB_DATA      = fullfile(paths.PIPELINE_ROOT, 'outputs', 'Lab_data');
paths.PSD_OUT       = fullfile(paths.PIPELINE_ROOT, 'outputs', 'psd_mat');
paths.POWER_OUT     = fullfile(paths.PIPELINE_ROOT, 'outputs', 'band_power');
paths.QC_PDF        = fullfile(paths.PIPELINE_ROOT, 'outputs', 'qc_pdf');
paths.QC_PNG        = fullfile(paths.PIPELINE_ROOT, 'outputs', 'qc_png');
paths.LOGS          = fullfile(paths.PIPELINE_ROOT, 'outputs', 'logs');
paths.DECISIONS     = fullfile(paths.PIPELINE_ROOT, 'outputs', 'decisions');
paths.MANIFEST_DIR  = fullfile(paths.PIPELINE_ROOT, 'outputs', 'manifests');
paths.METHODS_OUT   = fullfile(paths.PIPELINE_ROOT, 'outputs', 'methods');

% Channel location file (bundled with this repo).
paths.CHAN_LOC = fullfile(paths.PIPELINE_ROOT, 'chan_locs', 'Standard-10-20-Cap19.ced');

% EEGLAB installation. Auto-detected from the MATLAB path: add EEGLAB yourself
% (e.g. addpath('/path/to/eeglab'); eeglab nogui;) before running, or set
% EEGLAB_DIR explicitly in config_user.m.
eeglab_on_path = which('eeglab');
if ~isempty(eeglab_on_path)
    paths.EEGLAB_DIR = fileparts(eeglab_on_path);
else
    paths.EEGLAB_DIR = '';
end

% Optional downstream integration. Empty by default for new projects.
paths.DB             = '';
paths.PSY_ML_ROOT    = '';
paths.PSD_PYTHON_IN  = '';

% Per-project overrides. Copy/edit config_user.m for the senior's project.
if exist('config_user', 'file') == 2
    user_paths = config_user();
    names = fieldnames(user_paths);
    for i = 1:numel(names)
        paths.(names{i}) = user_paths.(names{i});
    end
end

% Auto-create output directories and the default raw_edf folder.
fields_to_mkdir = {'RAW_EDF_DIR','LAB_DATA','PSD_OUT','POWER_OUT','QC_PDF','QC_PNG','LOGS',...
                   'DECISIONS','MANIFEST_DIR','METHODS_OUT'};
for i = 1:numel(fields_to_mkdir)
    d = paths.(fields_to_mkdir{i});
    if ~isempty(d) && ~exist(d, 'dir')
        mkdir(d);
    end
end
end
