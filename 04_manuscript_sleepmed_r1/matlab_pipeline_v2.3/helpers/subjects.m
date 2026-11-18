function sids = subjects(opt)
% SUBJECTS  動態從 EDF 資料夾抓 S 開頭 ID(v2.2)
%
% USAGE
%   sids = subjects();           % 預設:全部(扣 EXCLUDED)
%   sids = subjects('all');      % 全部不扣
%   sids = subjects('excluded'); % 只回排除的
%   sids = subjects('count');    % 印出統計後仍回全部
%
% v2.2 ─ 2026/5/27

if nargin < 1; opt = 'default'; end

paths = config_paths();
pipe  = config_pipeline();

% 從 RAW_EDF_DIR 動態抓所有 S*.edf
files = dir(fullfile(paths.RAW_EDF_DIR, 'S*.edf'));
all_sids = sort({files.name});
all_sids = cellfun(@(x) x(1:end-4), all_sids, 'UniformOutput', false); % 去 .edf

excluded = pipe.excluded_subjects;

switch lower(opt)
    case {'all', 'raw'}
        sids = all_sids;
    case 'excluded'
        sids = intersect(all_sids, excluded);
    case 'count'
        n_total = numel(all_sids);
        n_excl  = numel(intersect(all_sids, excluded));
        fprintf('Subjects(來源 %s):\n', paths.RAW_EDF_DIR);
        fprintf('  總 .edf 檔: %d\n', n_total);
        fprintf('  排除清單命中: %d  (%s)\n', n_excl, strjoin(excluded, ', '));
        fprintf('  分析樣本: %d\n', n_total - n_excl);
        sids = setdiff(all_sids, excluded);
    otherwise   % 'default'
        sids = setdiff(all_sids, excluded);
end
end
