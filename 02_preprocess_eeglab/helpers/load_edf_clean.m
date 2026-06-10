function EEG = load_edf_clean(edf_path)
% LOAD_EDF_CLEAN  載入 EDF + 19-channel 標準化 + 套電極位置
%
% USAGE
%   EEG = load_edf_clean('/path/to/sub01.edf');
%
% 步驟:
%   1. pop_biosig 載入
%   2. 通道名標準化(去 -LE 後綴,case-fold)
%   3. 只保留 19 個標準通道
%   4. 重命名 T3/T4/T5/T6 → T7/T8/P7/P8(新 10-20)
%   5. 套用 Standard-10-20-Cap19.ced 之 3D 位置
%
% v2.2 ─ 2026/5/27

paths = config_paths();

EEG = pop_biosig(edf_path);
EEG = eeg_checkset(EEG);

% 通道名標準化對照表
std_map = { ...
    'FP1','Fp1'; 'FP2','Fp2'; ...
    'F7','F7';   'F3','F3';   'FZ','Fz';  'F4','F4';  'F8','F8'; ...
    'T3','T3';   'C3','C3';   'CZ','Cz';  'C4','C4';  'T4','T4'; ...
    'T5','T5';   'P3','P3';   'PZ','Pz';  'P4','P4';  'T6','T6'; ...
    'O1','O1';   'O2','O2' };

% 只保留 19 個標準通道
keep_idx = [];
for ci = 1:EEG.nbchan
    raw = EEG.chanlocs(ci).labels;
    clean = regexprep(upper(raw), '-.*$', '');  % 去 -LE / -A1 等
    match = find(strcmpi(std_map(:,1), clean));
    if ~isempty(match)
        EEG.chanlocs(ci).labels = std_map{match(1), 2};
        keep_idx(end+1) = ci; %#ok<AGROW>
    end
end
if length(keep_idx) < 19
    warning('load_edf_clean: 只找到 %d 個標準通道(預期 19)', length(keep_idx));
end
EEG = pop_select(EEG, 'channel', keep_idx);
EEG = eeg_checkset(EEG);

% T3/T4/T5/T6 → T7/T8/P7/P8(新 10-20 對應 .ced)
old2new = {'T3','T7'; 'T4','T8'; 'T5','P7'; 'T6','P8'};
for r = 1:size(old2new,1)
    idx = find(strcmp({EEG.chanlocs.labels}, old2new{r,1}));
    if ~isempty(idx)
        EEG.chanlocs(idx).labels = old2new{r,2};
    end
end

% 套用電極 3D 位置
locs_all = readlocs(paths.CHAN_LOC);
for ci = 1:EEG.nbchan
    m = find(strcmp({locs_all.labels}, EEG.chanlocs(ci).labels));
    if ~isempty(m)
        EEG.chanlocs(ci).X          = locs_all(m).X;
        EEG.chanlocs(ci).Y          = locs_all(m).Y;
        EEG.chanlocs(ci).Z          = locs_all(m).Z;
        EEG.chanlocs(ci).theta      = locs_all(m).theta;
        EEG.chanlocs(ci).radius     = locs_all(m).radius;
        EEG.chanlocs(ci).sph_theta  = locs_all(m).sph_theta;
        EEG.chanlocs(ci).sph_phi    = locs_all(m).sph_phi;
        EEG.chanlocs(ci).sph_radius = locs_all(m).sph_radius;
    end
end
EEG = eeg_checkset(EEG);
end
