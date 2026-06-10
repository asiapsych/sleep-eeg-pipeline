function run_pipeline_manual(raw_edf_dir, overwrite)
% RUN_PIPELINE_MANUAL  Preprocessing with human-in-the-loop confirmation.
%
% Same pipeline as RUN_PIPELINE_AUTO, but pauses for manual confirmation of
% the uncertain decisions instead of auto-confirming everything:
%   - review_ica_smart   : visually confirm mid-confidence ICLabel components
%                          (high/low confidence are auto-handled)
%   - review_epochs_smart: review epoch rejection for low-quality subjects
%
% This is the mode used for the manuscript: the machine runs the heavy
% steps, then a human confirms the ambiguous ICA components and epochs.
% Run this in an INTERACTIVE MATLAB session (it prompts for keypresses).
%
% USAGE
%   run_pipeline_manual()
%   run_pipeline_manual('/path/to/edf_folder')
%   run_pipeline_manual('/path/to/edf_folder', true)
%
% See also RUN_PIPELINE_AUTO (fully automated, no human confirmation).

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
    fprintf(fid, 'end\n');
    fclose(fid);
end

opts = struct();
opts.overwrite   = overwrite;
opts.skip_review = false;        % <-- human confirms ICA components & epochs
opts.start_from  = 'step1';
opts.stop_after  = 'step5';

subjects('count');
master_pipeline([], opts);
end
