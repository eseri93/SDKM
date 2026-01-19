%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double K-Means (DKM)
% Simultaneous partitioning (co-clustering) of rows and columns
%
% Reference: Maurizio Vichi (Oct 2012) - original MATLAB implementation
%
% Model: X ≈ U * Ym * V'  with hard partitions U (rows) and V (cols)
% Objective (equivalent forms):
%   - minimize ||X - U*Ym*V'||_F^2
%   - maximize ||U*Ym*V'||_F^2   (after centering/standardization)
%
% INPUTS
%   X        : (n x J) data matrix
%   K        : number of row clusters
%   Q        : number of column clusters
%   Rndstart : number of random starts (multistart)
%
% OUTPUTS
%   Vdkm   : (J x Q) hard membership matrix for variables (columns)
%   Udkm   : (n x K) hard membership matrix for units (rows)
%   Ymdkm  : (K x Q) centroid matrix
%   fdkm   : best explained variance (trace(B'B)/trace(Xs'Xs))
%   indkm  : iterations used in the best run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vdkm, Udkm, Ymdkm, fdkm, indkm] = DKM(X, K, Q, Rndstart)

%% ----------------------------
% 0) Basic checks and defaults
% ----------------------------
if nargin < 3
    error('DKM:TooFewInputs', 'Usage: DKM(X, K, Q, Rndstart).');
end
if nargin < 4 || isempty(Rndstart)
    Rndstart = 20;
end

if isempty(X) || ~isnumeric(X) || ndims(X) ~= 2
    error('DKM:InvalidX', 'X must be a non-empty numeric 2D matrix.');
end
[n, J] = size(X);

if ~isscalar(K) || K < 2 || K > n || K ~= round(K)
    error('DKM:InvalidK', 'K must be an integer in [2, n].');
end
if ~isscalar(Q) || Q < 2 || Q > J || Q ~= round(Q)
    error('DKM:InvalidQ', 'Q must be an integer in [2, J].');
end
if ~isscalar(Rndstart) || Rndstart < 1 || Rndstart ~= round(Rndstart)
    error('DKM:InvalidRndstart', 'Rndstart must be a positive integer.');
end

%% ----------------------------
% 1) Algorithm settings
% ----------------------------
maxIter = 100;        % maximum number of iterations per start
tol     = 1e-10;      % convergence tolerance on objective improvement

%% ----------------------------
% 2) Centering and standardization
% ----------------------------
% Row-centering matrix
Jm = eye(n) - (1/n) * ones(n);

% Column variances after row-centering
S  = (1/n) * (X' * Jm * X);
v  = diag(S);

% Safeguard: avoid division by zero if some columns have zero variance
epsVar = 1e-12;
v(v < epsVar) = epsVar;

% Standardized data (row-centered, column-scaled)
Xs = Jm * X * diag(1 ./ sqrt(v));

% Total deviance (Frobenius norm squared)
st = sum(sum(Xs.^2));

%% ----------------------------
% 3) Multistart loop
% ----------------------------
fdkm   = -Inf;
Vdkm   = [];
Udkm   = [];
Ymdkm  = [];
indkm  = 0;

for loop = 1:Rndstart

    % Random hard partitions (non-empty by construction)
    V = randPU(J, Q);
    U = randPU(n, K);

    su = sum(U, 1);
    sv = sum(V, 1);

    % Initial centroid matrix
    Ym = diag(1./su) * (U' * Xs) * V * diag(1./sv);

    % Initial fit (explained variance)
    B  = U * Ym * V';
    f0 = trace(B' * B) / st;

    fdif = Inf;
    it   = 0;

    %% Iteration phase
    while (fdif > tol) && (it < maxIter)
        it = it + 1;

        % ----------------------------
        % (A) Update U given (Ym, V)
        % ----------------------------
        U = zeros(n, K);
        Ymv = Ym * V';  % (K x J) row-cluster prototypes in feature space

        for i = 1:n
            % Assign row i to the closest row-prototype (Euclidean distance)
            mindif = sum((Xs(i,:) - Ymv(1,:)).^2);
            posmin = 1;
            for k = 2:K
                dif = sum((Xs(i,:) - Ymv(k,:)).^2);
                if dif < mindif
                    mindif = dif;
                    posmin = k;
                end
            end
            U(i, posmin) = 1;
        end

        % Ensure non-empty row clusters (rebalance if needed)
        U = balanceCluster(U);

        su = sum(U, 1);

        % ----------------------------
        % (B) Update Ym given (U, V)
        % ----------------------------
        Ym = diag(1./su) * (U' * Xs) * V * diag(1./sv);

        % ----------------------------
        % (C) Update V given (Ym, U)
        % ----------------------------
        V = zeros(J, Q);
        Ymu = U * Ym;   % (n x Q) column-cluster prototypes in observation space

        for j = 1:J
            % Assign column j to the closest column-prototype (Euclidean distance)
            mindif = sum((Xs(:,j) - Ymu(:,1)).^2);
            posmin = 1;
            for q = 2:Q
                dif = sum((Xs(:,j) - Ymu(:,q)).^2);
                if dif < mindif
                    mindif = dif;
                    posmin = q;
                end
            end
            V(j, posmin) = 1;
        end

        % Ensure non-empty column clusters
        V = balanceCluster(V);

        sv = sum(V, 1);

        % ----------------------------
        % (D) Update Ym again given updated (U, V)
        % ----------------------------
        Ym = diag(1./su) * (U' * Xs) * V * diag(1./sv);

        % ----------------------------
        % (E) Objective and stopping rule
        % ----------------------------
        B  = U * Ym * V';
        f  = trace(B' * B) / st;

        fdif = abs(f - f0);
        f0   = f;
    end

    fprintf('DKM: loop=%d | fit=%.6g | iter=%d | fdif=%.3g\n', loop, f0, it, fdif);

    % Keep best run
    if f0 > fdkm
        fdkm  = f0;
        Udkm  = U;
        Vdkm  = V;
        Ymdkm = Ym;
        indkm = it;
    end
end

%% ----------------------------
% 4) Sort clusters by descending cardinality (and reorder Ym consistently!)
% ----------------------------

% Sort column clusters (Q) by size
[~, ordQ] = sort(sum(Vdkm,1), 'descend');
Vdkm  = Vdkm(:, ordQ);
Ymdkm = Ymdkm(:, ordQ);

% Sort row clusters (K) by size
[~, ordK] = sort(sum(Udkm,1), 'descend');
Udkm  = Udkm(:, ordK);
Ymdkm = Ymdkm(ordK, :);

fprintf('DKM (Final): fit=%.6g | best-iter=%d | best-loop (saved)\n', fdkm, indkm);

end % <-- end main function


%% ============================================================
% Local utilities
%% ============================================================

function U = randPU(n, c)
%RANDPU Generate a random hard partition of n items into c non-empty classes.
U = zeros(n, c);
U(1:c, :) = eye(c);
U(c+1:n, 1) = 1;
for i = c+1:n
    U(i, 1:c) = U(i, randperm(c));
end
U = U(randperm(n), :);
end


function A = balanceCluster(A)
%BALANCECLUSTER Ensure that no cluster is empty by splitting the largest one.
sv = sum(A, 1);
while any(sv == 0)
    [~, p1] = min(sv); % empty (or smallest) cluster
    [~, p2] = max(sv); % largest cluster
    ind = find(A(:, p2));
    ind = ind(1:floor(sv(p2)/2));
    A(ind, p1) = 1;
    A(ind, p2) = 0;
    sv = sum(A, 1);
end
end
