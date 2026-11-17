function paths = config_paths()
% CONFIG_PATHS  集中管理所有檔案路徑（v2.3，公開版）
%
% 研究時這個檔寫的是本機的絕對路徑。公開版改成由環境變數指定；沒有設定時用本資料夾底下的預設位置。
%   EEG_RAW_EDF_DIR   原始 EDF 所在的資料夾（檔名為 <受試者編號>.edf）
%   EEG_RAW_NG_DIR    NeuroGuide 原始檔所在的資料夾（只有轉檔時用得到）
%   EEGLAB_DIR        EEGLAB 的安裝位置（未設定時用 MATLAB 路徑上找到的 eeglab）
% 其餘輸出都寫在本資料夾的 outputs/ 底下。

paths = struct();
root = fileparts(fileparts(mfilename('fullpath')));      % matlab_pipeline_v2.3/

% ── 輸入 ──
paths.RAW_EDF_DIR = env_or('EEG_RAW_EDF_DIR', fullfile(root, 'data', 'raw_edf'));
paths.RAW_NG_DIR  = env_or('EEG_RAW_NG_DIR',  fullfile(root, 'data', 'raw_ng'));
paths.NG_MANIFEST = fullfile(paths.RAW_NG_DIR, '_manifest.csv');

% ── 工作目錄（EEGLAB 處理的中間結果）──
paths.PIPELINE_ROOT = root;
paths.LAB_DATA      = fullfile(paths.PIPELINE_ROOT, 'outputs', 'Lab_data');
paths.PSD_OUT       = fullfile(paths.PIPELINE_ROOT, 'outputs', 'psd_mat');
paths.QC_PDF        = fullfile(paths.PIPELINE_ROOT, 'outputs', 'qc_pdf');
paths.QC_PNG        = fullfile(paths.PIPELINE_ROOT, 'outputs', 'qc_png');
paths.LOGS          = fullfile(paths.PIPELINE_ROOT, 'outputs', 'logs');
paths.DECISIONS     = fullfile(paths.PIPELINE_ROOT, 'outputs', 'decisions');
paths.MANIFEST_DIR  = fullfile(paths.PIPELINE_ROOT, 'outputs', 'manifests');
paths.METHODS_OUT   = fullfile(paths.PIPELINE_ROOT, 'outputs', 'methods');

% ── 通道位置檔 ──
paths.CHAN_LOC = fullfile(root, 'chan_locs', 'Standard-10-20-Cap19.ced');

% ── EEGLAB 安裝位置 ──
paths.EEGLAB_DIR = env_or('EEGLAB_DIR', fileparts(which('eeglab')));

% ── 自動建立所有 output 目錄 ──
fields_to_mkdir = {'LAB_DATA','PSD_OUT','QC_PDF','QC_PNG','LOGS',...
                   'DECISIONS','MANIFEST_DIR','METHODS_OUT'};
for i = 1:numel(fields_to_mkdir)
    d = paths.(fields_to_mkdir{i});
    if ~exist(d, 'dir')
        mkdir(d);
    end
end
end


function v = env_or(name, default)
v = getenv(name);
if isempty(v)
    v = default;
end
end
