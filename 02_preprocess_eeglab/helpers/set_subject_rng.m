function seed = set_subject_rng(sid, global_seed)
% SET_SUBJECT_RNG  以 subject ID 為基礎產生 deterministic seed
%
% USAGE
%   seed = set_subject_rng('sub01');             % 預設 global_seed=42
%   seed = set_subject_rng('sub01', 42);         % 顯式 global_seed
%
% 設計:
%   - hash subject ID → 32-bit unsigned int → 跟 global_seed 結合
%   - 同人 hash 永遠相同 → reproducibility
%   - 跨人不同 hash → 隨機操作獨立
%
% 範例輸出(deterministic):
%   set_subject_rng('sub01') → 某固定 uint32
%
% 也會呼叫 rng() 設定當前 random state。
%
% v2.2 ─ 2026/5/27

if nargin < 2
    pipe = config_pipeline();
    global_seed = pipe.global_seed;
end

% 簡單但 deterministic 的 hash:用 char codes 加總 + global_seed
sid_str = char(sid);
hash = uint32(global_seed);
for i = 1:length(sid_str)
    hash = uint32(mod(double(hash) * 31 + double(sid_str(i)), 2^32));
end
seed = hash;

% 設定當前 random state
rng(double(seed), 'twister');
end
