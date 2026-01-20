function [Vdskm, Udskm, Ymdskm, fdskm, indskm] = SDKM(X, K, Q, varargin)
%SDKM  Spherical Double k-Means (SDKM) co-clustering algorithm.
%
% Reference: Bombelli, I., Iezzi, D.F., Seri, E., & Vichi, M. 
% "Spherical Double K-Means: a co-clustering approach for textual data analysis".
%
%   [V,U,Ym,f,iterBest] = SDKM(X,K,Q,...) simultaneously partitions:
%     - rows (objects/documents) into K clusters via U (n x K)
%     - columns (variables/terms) into Q clusters via V (J x Q)
%
%   INPUTS
%     X  : (n x J) data matrix (objects x variables)
%     K  : number of object clusters
%     Q  : number of variable clusters
%
%   NAME-VALUE OPTIONS
%     'Stats'    : 'on' (default) or 'off' -> print iteration summaries
%     'Stand'    : 'on' (default), 'off', 'Mahalanobis', 'DSt' (double-std)
%                 (alias 'Dst' is also accepted)
%     'Rndst'    : number of random starts (default 20)
%     'MaxIter'  : maximum iterations per start (default 100)
%     'ConvToll' : convergence tolerance on |f - f0| (default 1e-9)

%% ----------------------------
% 0) Basic input checks
% ----------------------------
if nargin < 3
    error('SDKM:TooFewInputs', 'Usage: SDKM(X, K, Q, ...)');
end

if isempty(X) || ~isnumeric(X)
    error('SDKM:InvalidX', 'X must be a non-empty numeric matrix.');
end
if min(size(X)) == 1
    error('SDKM:NotEnoughData', 'X must be a 2D matrix with n>1 and J>1.');
end

[n, J] = size(X);

if isempty(K) || ~isnumeric(K) || ~isscalar(K) || K < 1 || K > n || K ~= floor(K)
    error('SDKM:InvalidK', 'K must be an integer in [1, n].');
end
if isempty(Q) || ~isnumeric(Q) || ~isscalar(Q) || Q < 1 || Q > J || Q ~= floor(Q)
    error('SDKM:InvalidQ', 'Q must be an integer in [1, J].');
end

%% ----------------------------
% 1) Parse optional parameters
% ----------------------------
p = inputParser;
p.FunctionName = 'SDKM';

addParameter(p, 'Stats',    'on',  @(s) ischar(s) || isstring(s));
addParameter(p, 'Stand',    'on',  @(s) ischar(s) || isstring(s));
addParameter(p, 'Rndst',    20,    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1000);
addParameter(p, 'MaxIter',  100,   @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1000);
addParameter(p, 'ConvToll', 1e-9,  @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 0.1);

parse(p, varargin{:});
Stats    = lower(string(p.Results.Stats));
Stand    = string(p.Results.Stand);
Rndst    = p.Results.Rndst;
MaxIter  = p.Results.MaxIter;
ConvToll = p.Results.ConvToll;

if ~(Stats == "on" || Stats == "off")
    error('SDKM:InvalidStats', 'Stats must be ''on'' or ''off''.');
end

% Accept a couple of variants for double standardization naming
Stand = strrep(Stand, "Dst", "DSt");
Stand = strrep(Stand, "dst", "DSt");

%% ----------------------------
% 2) Standardize / preprocess X
% ----------------------------
StandLower = lower(Stand);
switch StandLower
    case "off"
        % Center columns only
        Xs = X - ones(n,1) * mean(X, 1);
    case "on"
        % Z-score columns (population normalization)
        Xs = zscore(X, 1);
    case "mahalanobis"
        % Whitening-like transform (may be unstable if cov is singular)
        Jc = eye(n) - (1/n) * ones(n);
        Sx = cov(X, 1);
        Xs = Jc * X * (Sx^(-0.5));
    case "dst"
        % Double standardization (iterative row+column z-scoring)
        [Xs, ~, ~] = DStand(X);
        Xs = Xs ./ sqrt(n);
    otherwise
        error('SDKM:InvalidStand', ...
            'Stand must be ''on'', ''off'', ''Mahalanobis'', or ''DSt''.');
end

%% ----------------------------
% 3) Multi-start optimization
% ----------------------------
UC = eye(K);
VC = eye(Q);
KK = 1:K;
QQ = 1:Q;

fdskm    = -Inf;
indskm   = 0;
loopBest = 1;
fdifBest = NaN;

Xs_norm = norm(Xs, 'fro');
epsDen  = 1e-12;

for loop = 1:Rndst

    % Random non-empty hard partitions
    U = randPU(n, K);   % (n x K)
    V = randPU(J, Q);   % (J x Q)

    % Initial model matrices
    [Ym, Xt, f0] = updateYmAndObjective(Xs, U, V, epsDen);

    % Iterate updates
    for it = 1:MaxIter

        % ---- Update U (object clusters) given current structure ----
        % Use cosine-like assignment against current aggregated profile
        dD  = sqrt(diag((U' * Xs) * (Xs' * U)));
        dD  = max(dD, epsDen);
        Dm1 = diag(1 ./ dD);

        Ymv = Dm1 * U' * Xs;  % (K x J)

        for i = 1:n
            posmax = KK(U(i,:) == 1);

            den0   = sqrt( (Xs(i,:) * Xs(i,:)') * (Ymv(posmax,:) * Ymv(posmax,:)' ) );
            den0   = max(den0, epsDen);
            maxdif = (Xs(i,:) * Ymv(posmax,:)' ) / den0;

            if sum(U(:,posmax)) > 1
                for k = find(KK ~= posmax)
                    U(i,:) = UC(k,:);
                    denk   = sqrt( (Xs(i,:) * Xs(i,:)') * (Ymv(k,:) * Ymv(k,:)' ) );
                    denk   = max(denk, epsDen);
                    dif    = (Xs(i,:) * Ymv(k,:)' ) / denk;

                    if dif > maxdif
                        maxdif = dif;
                        posmax = k;
                    end
                end
                U(i,:) = UC(posmax,:);
            end
        end
        U = balanceCluster(U);

        % ---- Update Ym and objective ----
        [Ym, Xt, f] = updateYmAndObjective(Xs, U, V, epsDen);

        % ---- Update V (variable clusters) ----
        dF  = sqrt(diag((V' * Xs') * (Xs * V)));
        dF  = max(dF, epsDen);
        Fm1 = diag(1 ./ dF);

        Ymu = Xs * V * Fm1;  % (n x Q)

        for j = 1:J
            posmax = QQ(V(j,:) == 1);

            den0   = sqrt( (Xs(:,j)' * Xs(:,j)) * (Ymu(:,posmax)' * Ymu(:,posmax)) );
            den0   = max(den0, epsDen);
            maxdif = (Xs(:,j)' * Ymu(:,posmax)) / den0;

            if sum(V(:,posmax)) > 1
                for q = find(QQ ~= posmax)
                    V(j,:) = VC(q,:);
                    denq   = sqrt( (Xs(:,j)' * Xs(:,j)) * (Ymu(:,q)' * Ymu(:,q)) );
                    denq   = max(denq, epsDen);
                    dif    = (Xs(:,j)' * Ymu(:,q)) / denq;

                    if dif > maxdif
                        maxdif = dif;
                        posmax = q;
                    end
                end
                V(j,:) = VC(posmax,:);
            end
        end
        V = balanceCluster(V);

        % ---- Update Ym and objective again ----
        [Ym, Xt, f] = updateYmAndObjective(Xs, U, V, epsDen);

        fdif = abs(f - f0);
        if fdif > ConvToll
            f0 = f;
        else
            break;
        end
    end

    if Stats == "on"
        fprintf('SDKM: start=%d, explainedVar=%g, iter=%d, fdif=%g\n', loop, f, it, fdif);
    end

    % Keep best solution over random starts
    if loop == 1 || f > fdskm
        Vdskm    = V;
        Udskm    = U;
        Ymdskm   = Ym;
        fdskm    = f;
        indskm   = it;
        loopBest = loop;
        fdifBest = fdif;
    end
end

%% ----------------------------
% 4) Post-processing: sort clusters by size
% ----------------------------
[~, icV] = sort(diag(Vdskm' * Vdskm), 'descend');
Vdskm = Vdskm(:, icV);

[~, icU] = sort(diag(Udskm' * Udskm), 'descend');
Udskm = Udskm(:, icU);

if Stats == "on"
    fprintf('SDKM (final): explainedVar=%g, bestStart=%d, iter=%d, fdif=%g\n', ...
        fdskm, loopBest, indskm, fdifBest);
end

end % ===== end main function =====


%% ========================================================================
% Helper: update Ym, reconstruct Xt, compute objective f
% ========================================================================
function [Ym, Xt, f] = updateYmAndObjective(Xs, U, V, epsDen)
    % Cluster sizes
    dU = diag(U' * U);
    dV = diag(V' * V);
    dU = max(dU, epsDen);
    dV = max(dV, epsDen);

    Um1 = diag(1 ./ dU);
    Vm1 = diag(1 ./ dV);

    % Core association matrix (K x Q), then row-normalize (spherical)
    Ym = Um1 * (U' * Xs) * V * Vm1;
    Ym = rowNormalize(Ym, epsDen);

    % Reconstruction in data space
    Xt = U * Ym * V';

    % Explained-variance-like objective (cosine between Xs and Xt in Fro space)
    num = sum(sum(Xs .* Xt));
    den = max(norm(Xs, 'fro') * norm(Xt, 'fro'), epsDen);
    f   = num / den;
end


%% ========================================================================
% Helper: row normalization (avoid dependency on normr)
% ========================================================================
function A = rowNormalize(A, epsDen)
    rn = sqrt(sum(A.^2, 2));
    rn = max(rn, epsDen);
    A  = A ./ rn;
end


%% ========================================================================
% Helper: random non-empty hard partition matrix (n x c)
% ========================================================================
function U = randPU(n, c)
    % Generates a random partition of n items into c non-empty classes.
    U = zeros(n, c);
    U(1:c, :) = eye(c);
    U(c+1:n, 1) = 1;

    for i = c+1:n
        U(i, 1:c) = U(i, randperm(c));
    end

    U = U(randperm(n), :);
end


%% ========================================================================
% Helper: double standardization (iterative column+row z-scoring)
% ========================================================================
function [Xs, it, f] = DStand(X)
    epsTol = 1e-14;
    Xso = X;

    for it = 1:30000
        X1 = zscore(Xso, 1);      % standardize columns
        Xs = (zscore(X1', 1))';   % standardize rows
        f  = norm(Xs - Xso);

        if f < epsTol
            break;
        end
        Xso = Xs;
    end
end


%% ========================================================================
% Helper: ensure no empty clusters by rebalancing
% ========================================================================
function A = balanceCluster(A)
    sv = sum(A, 1);
    while any(sv == 0)
        [~, p1] = min(sv);
        [~, p2] = max(sv);

        ind = find(A(:, p2));
        ind = ind(1:floor(sv(p2) / 2));

        A(ind, p1) = 1;
        A(ind, p2) = 0;
        sv = sum(A, 1);
    end
end

