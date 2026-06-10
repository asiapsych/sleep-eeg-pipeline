function run_pipeline_auto(raw_edf_dir, overwrite)
% RUN_PIPELINE_AUTO  One-command fully automated preprocessing and power export.
%
% USAGE
%   run_pipeline_auto()
%   run_pipeline_auto('/path/to/new/edf_folder')
%   run_pipeline_auto('/path/to/new/edf_folder', true)

if nargin < 1
    raw_edf_dir = '';
end
if nargin < 2
    overwrite = false;
end

root = fileparts(mfilename('fullpath'));
addpath(genpath(root));

if ~isempty(raw_edf_dir)
    cfg_path = fullfile(root, 'config', 'config_user.m');
    fid = fopen(cfg_path, 'w');
    fprintf(fid, 'function paths = config_user()\n');
    fprintf(fid, 'paths = struct();\n');
    fprintf(fid, 'paths.RAW_EDF_DIR = ''%s'';\n', strrep(raw_edf_dir, '''', ''''''));
    fprintf(fid, 'paths.RAW_FILE_PATTERN = ''*.edf'';\n');
    fprintf(fid, 'paths.DATASET_LABEL = ''senior_collection'';\n');
    fprintf(fid, 'end\n');
    fclose(fid);
end

opts = struct();
opts.overwrite = overwrite;
opts.skip_review = true;
opts.start_from = 'step1';
opts.stop_after = 'step5';

subjects('count');
master_pipeline([], opts);
end
