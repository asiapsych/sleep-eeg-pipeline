function bands = config_bands()
% CONFIG_BANDS  頻段定義(v2.2)─ 對齊 task_f_recompute_relpower.py
%
% USAGE
%   b = config_bands();
%   b.primary.Beta   % [13 30]
%   b.derived.Beta1  % [13 18]
%
% PRIMARY:5 個非重疊連續頻段(進 relative power 分母)
% DERIVED:sub-band(僅算 子/total,不進分母)
%
% v2.2 ─ 2026/5/27

bands = struct();

% ── PRIMARY(進分母,非重疊連續)──
bands.primary.Delta = [1.0  4.0];
bands.primary.Theta = [4.0  8.0];
bands.primary.Alpha = [8.0  13.0];
bands.primary.Beta  = [13.0 30.0];
bands.primary.Gamma = [30.0 45.0];

% ── DERIVED(只算 子/total,不進分母)──
bands.derived.Alpha1   = [8.0  10.0];
bands.derived.Alpha2   = [10.0 13.0];
bands.derived.Beta1    = [13.0 18.0];  % 本研究定義(paper Methods 註明)
bands.derived.Beta2    = [18.0 24.0];
bands.derived.Beta3    = [24.0 30.0];
bands.derived.HighBeta = [30.0 40.0];

% ── 全部合併(供 PSD 圖、QC 用)──
bands.all = struct();
primary_names = fieldnames(bands.primary);
for i = 1:numel(primary_names)
    bands.all.(primary_names{i}) = bands.primary.(primary_names{i});
end
derived_names = fieldnames(bands.derived);
for i = 1:numel(derived_names)
    bands.all.(derived_names{i}) = bands.derived.(derived_names{i});
end

% ── 報告用:文獻對應(paper Methods 註明)──
bands.notes = struct();
bands.notes.beta_subband_caveat = ...
    ['本研究 Beta 子頻段定義為 Beta1=13-18, Beta2=18-24, Beta3=24-30 Hz。' ...
     '文獻中 Beta sub-band 定義不一,常見之「Beta2 = 24-32 Hz」hyperarousal marker ' ...
     '在本研究中對應 Beta3 區段(24-30 Hz)。'];
end
