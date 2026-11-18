function T = surface_laplacian_perrin(xyz, stiffness, lambda, n_terms)
% SURFACE_LAPLACIAN_PERRIN  球面樣條表面拉普拉斯（surface Laplacian）的轉換矩陣
%
% USAGE
%   T   = surface_laplacian_perrin(xyz);          % xyz: [n_ch × 3] 電極座標
%   lap = T * data;                               % data: [n_ch × n_samples]
%
% 做法：Perrin et al. (1989, Electroencephalogr Clin Neurophysiol 72:184–187；1990 更正) 的球面樣條。
% 參數的預設值與 MNE-Python 的 compute_current_source_density 相同：
%   stiffness = 4、lambda = 1e-5、Legendre 多項式取 50 項。
% 電極座標先投影到單位球面。輸出的單位是任意的（未除以頭半徑的平方）；
% 後續只使用對數頻譜上的相對量（峰值功率、去背景功率、相對功率），不受整體尺度影響。
%
% 表面拉普拉斯與參考電極無關，會壓低空間上分布廣的成分（容積傳導、遠處的肌電），
% 保留電極下方局部的活動。
%
% v2.3 ─ 2026/11/17

if nargin < 2 || isempty(stiffness); stiffness = 4; end
if nargin < 3 || isempty(lambda);    lambda = 1e-5; end
if nargin < 4 || isempty(n_terms);   n_terms = 50; end

n = size(xyz, 1);
xyz = xyz ./ sqrt(sum(xyz .^ 2, 2));                  % 投影到單位球面
cosang = max(min(xyz * xyz', 1), -1);                 % 兩兩電極夾角的餘弦

% Legendre 多項式（遞迴式）：P0 = 1、P1 = x、(k+1) P_{k+1} = (2k+1) x P_k − k P_{k−1}
G = zeros(n); H = zeros(n);
Pkm1 = ones(n); Pk = cosang;
for k = 1:n_terms
    if k > 1
        Pkp1 = ((2 * k - 1) .* cosang .* Pk - (k - 1) .* Pkm1) ./ k;
        Pkm1 = Pk; Pk = Pkp1;
    end
    G = G + (2 * k + 1) / (k * (k + 1)) ^ stiffness       .* Pk;
    H = H + (2 * k + 1) / (k * (k + 1)) ^ (stiffness - 1) .* Pk;
end
G = G / (4 * pi); H = H / (4 * pi);

% 樣條係數（含常數項的限制式），再乘上 H 得到拉普拉斯
Gi  = inv(G + lambda * eye(n));
TC  = sum(Gi, 1);                                     % 1 × n
sgi = sum(Gi(:));
Z   = eye(n) - 1 / n;                                 % 先減去各時間點的電極平均
Cp  = Gi * Z;
C   = Cp - (TC' * sum(Cp, 1)) / sgi;
T   = H' * C;
end
