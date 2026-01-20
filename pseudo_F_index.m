function pF_dk = pseudo_F_index(X, U, V)
%PSEUDO_F_INDEX  Pseudo-F index for SDKM co-clustering solutions.
%
% Reference: Bombelli, I., Iezzi, D.F., Seri, E., & Vichi, M. 
% "Spherical Double K-Means: a co-clustering approach for textual data analysis".
%
%   pF_dk = pseudo_F_index(X, U, V) computes the pseudo-F index based on
%   Frobenius-norm between/within deviance under the co-clustering
%   projection defined by U (rows) and V (columns).
%
%   INPUTS
%     X : (n x J) data matrix
%     U : (n x K) hard membership matrix for units (rows)
%     V : (J x Q) hard membership matrix for variables (columns)
%
%   OUTPUT
%     pF_dk : pseudo-F value (larger suggests better separation)
%
%   Notes:
%   - Assumes U and V represent non-empty hard partitions (one 1 per row).
%   - Uses projection matrices H_U = U (U'U)^{-1} U' and H_V similarly.

%% ----------------------------
% 0) Basic checks
% ----------------------------
if nargin < 3
    error('pseudo_F_index:TooFewInputs', 'Usage: pseudo_F_index(X, U, V).');
end

if isempty(X) || ~isnumeric(X) || ndims(X) ~= 2
    error('pseudo_F_index:InvalidX', 'X must be a non-empty numeric 2D matrix.');
end

[n, J] = size(X);

if isempty(U) || ~isnumeric(U) || size(U,1) ~= n
    error('pseudo_F_index:InvalidU', 'U must be numeric with size(U,1)=size(X,1).');
end
if isempty(V) || ~isnumeric(V) || size(V,1) ~= J
    error('pseudo_F_index:InvalidV', 'V must be numeric with size(V,1)=size(X,2).');
end

K = size(U,2);
Q = size(V,2);

if K < 1 || Q < 1
    error('pseudo_F_index:InvalidKQ', 'U and V must have at least one cluster (K>=1, Q>=1).');
end

% Check for empty clusters (columns with zero members)
su = sum(U, 1);
sv = sum(V, 1);
if any(su == 0) || any(sv == 0)
    error('pseudo_F_index:EmptyCluster', 'U and/or V contains empty clusters.');
end

%% ----------------------------
% 1) Projection matrices
% ----------------------------
% For hard partitions, U'*U and V'*V are diagonal (cluster sizes),
% but we keep a generic implementation with a small safeguard.
epsDen = 1e-12;

UtU = U' * U;
VtV = V' * V;

% Safeguard against numerical issues
UtU = UtU + epsDen * eye(size(UtU));
VtV = VtV + epsDen * eye(size(VtV));

H_U = U / UtU * U';   % (n x n)
H_V = V / VtV * V';   % (J x J)

%% ----------------------------
% 2) Grand mean matrix
% ----------------------------
one_n = ones(n, 1);
one_J = ones(J, 1);

% Equivalent to filling an (n x J) matrix with the global mean of X
grand_mean_matrix = (one_n * one_n') * X * (one_J * one_J') / (n * J);

%% ----------------------------
% 3) Between/within deviance
% ----------------------------
X_hat = H_U * X * H_V;

SS_between = norm(X_hat - grand_mean_matrix, 'fro')^2;
SS_within  = norm(X - X_hat, 'fro')^2;

%% ----------------------------
% 4) Degrees of freedom and pseudo-F
% ----------------------------
df_between = (K * Q) - 1;
df_within  = (n * J) - (K * Q);

if df_between <= 0 || df_within <= 0
    error('pseudo_F_index:InvalidDF', ...
        'Invalid degrees of freedom: check that K*Q is in [2, n*J-1].');
end

if SS_within <= epsDen
    % Perfect reconstruction (or numerical zero within deviance)
    pF_dk = Inf;
else
    pF_dk = (SS_between / df_between) / (SS_within / df_within);
end

end
