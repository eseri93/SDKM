function [X, Ym, U_true, V_true] = GenDataSDKM(n, K, J, Q, er1, er2, varargin)
%GENDATASDKM  Synthetic data generator for Spherical Double K-Means (SDKM).
%
%   [X, Ym, U_true, V_true] = GenDataSDKM(n, K, J, Q, er1, er2)
%   generates a synthetic (n x J) data matrix X from a co-clustering model:
%
%       X = U_true * Ym * V_true' + E
%
%   where:
%     - U_true is an (n x K) hard partition of rows (objects/documents),
%     - V_true is a (J x Q) hard partition of columns (variables/terms),
%     - Ym is a (K x Q) matrix of unit-norm row centroids in spherical space,
%     - E is Gaussian noise with std = er2.
%
%   INPUTS
%     n   : number of objects (rows)
%     K   : number of row clusters
%     J   : number of variables (columns)
%     Q   : number of column clusters (must satisfy 2 <= Q <= min(K-1, J))
%     er1 : centroid perturbation level (heterogeneity of centroids)
%     er2 : noise level (heterogeneity within clusters)
%
%   OPTIONAL NAME-VALUE PAIRS
%     'Plot' : 'on'/'off' (default 'off') - show imagesc(X)
%     'Seed' : scalar integer (default []) - RNG seed for reproducibility
%
%   OUTPUTS
%     X      : (n x J) generated data matrix
%     Ym     : (K x Q) row-normalized centroid matrix
%     U_true : (n x K) true hard partition (rows)
%     V_true : (J x Q) true hard partition (columns)

%% ----------------------------
% 0) Optional arguments
% ----------------------------
p = inputParser;
p.FunctionName = 'GenDataSDKM';
addParameter(p, 'Plot', 'off', @(s) ischar(s) || isstring(s));
addParameter(p, 'Seed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
parse(p, varargin{:});

plotFlag = lower(string(p.Results.Plot));
seedVal  = p.Results.Seed;

if ~isempty(seedVal)
    rng(seedVal);
end

%% ----------------------------
% 1) Input checks
% ----------------------------
if ~(isscalar(n) && n == floor(n) && n > 0), error('n must be a positive integer.'); end
if ~(isscalar(K) && K == floor(K) && K > 1), error('K must be an integer >= 2.'); end
if ~(isscalar(J) && J == floor(J) && J > 0), error('J must be a positive integer.'); end
if ~(isscalar(Q) && Q == floor(Q) && Q > 1), error('Q must be an integer >= 2.'); end

if Q > min(K-1, J)
    error('Q must satisfy Q <= min(K-1, J). Given K=%d, J=%d, Q=%d.', K, J, Q);
end

if ~(isscalar(er1) && isnumeric(er1) && er1 >= 0), error('er1 must be a non-negative scalar.'); end
if ~(isscalar(er2) && isnumeric(er2) && er2 >= 0), error('er2 must be a non-negative scalar.'); end

%% ----------------------------
% 2) Generate true partitions
% ----------------------------
U_true = randPU(n, K);    % (n x K)
V_true = randPU(J, Q);    % (J x Q)

% Sort rows/columns by cluster label for nicer block-structure visualization
[~, idxV] = sort(V_true * (1:Q)');   % cluster label per variable
V_true = V_true(idxV, :);

[~, idxU] = sort(U_true * (1:K)');   % cluster label per object
U_true = U_true(idxU, :);

%% ----------------------------
% 3) Build a "roughly equidistant" centroid configuration in K x Q
% ----------------------------
% Construct a simple distance matrix among K row-clusters:
D  = 10 * (ones(K) - eye(K));        % zero on diagonal, constant off-diagonal
Jc = eye(K) - (1 / K) * ones(K);     % centering matrix
B  = -0.5 * Jc * D * Jc;             % doubly centered matrix (Gram-like)

% Extract Q leading components (B is symmetric)
[EigVec, EigVal] = eigs(B, Q, 'la');

% Defensive: keep only non-negative eigenvalues (numerical safeguard)
lam = diag(EigVal);
lam(lam < 0) = 0;
Ym = EigVec * diag(sqrt(lam));

% Add centroid perturbation (controls centroid heterogeneity) and normalize rows
Ym = Ym + randn(K, Q) * er1;
Ym = rowNormalize(Ym);

%% ----------------------------
% 4) Generate data matrix
% ----------------------------
X = U_true * Ym * V_true' + randn(n, J) * er2;

%% ----------------------------
% 5) Optional visualization
% ----------------------------
if plotFlag == "on"
    imagesc(X);
    colormap(flipud(hot));
    colorbar;
    title('Generated SDKM data matrix X');
end

end

%% ========================================================================
% Local helpers
% ========================================================================

function U = randPU(n, c)
%RANDPU  Random hard partition with non-empty clusters.
%
%   U = randPU(n, c) returns an (n x c) matrix with exactly one 1 per row
%   and no empty columns.

if c > n
    error('randPU:InvalidSizes', 'Number of clusters c cannot exceed n.');
end

U = zeros(n, c);
U(1:c, :) = eye(c);          % ensure non-empty clusters
U(c+1:n, 1) = 1;             % fill remaining rows (temporarily)

for i = c+1:n
    U(i, 1:c) = U(i, randperm(c));
end

U = U(randperm(n), :);       % shuffle rows
end

function A = rowNormalize(A)
%ROWNORMALIZE  Normalize each row to have unit L2 norm (safe for zero rows).
rowNorms = sqrt(sum(A.^2, 2));
rowNorms(rowNorms == 0) = 1;
A = A ./ rowNorms;
end






    


