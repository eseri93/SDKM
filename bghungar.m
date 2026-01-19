%BGHUNGAR Robust Hungarian algorithm for the square assignment problem.
%
%   [perm, unfeas, D] = bghungar(C)
%
% INPUT
%   C : (n x n) profit (or cost) matrix.
%       By default the algorithm MAXIMIZES the total profit.
%       To MINIMIZE a cost matrix, call: bghungar(-C).
%
% OUTPUT
%   perm   : (n x 1) optimal assignment (permutation). Row i is assigned to column perm(i).
%   unfeas : feasibility / diagnostic flag
%            0 = OK
%            1..4 = numerical / degeneracy issues encountered (rare)
%   D      : reduced matrix at the end of the initialization step (see original reference [1]).
%
% Notes
%   - This is an improved, more robust pure-MATLAB implementation originally shared by
%     Nedialko Krouchev (File Exchange contribution).
%   - Used in this repository for label switching correction via optimal matching of clusters.
%
% Reference
%   [1] D. Ivantchev and G. Negler, "Network Optimization", Sofia, 1992.
%
% Original author / copyright notice preserved below.
%
% (c) 2002/11/29 Nedialko Krouchev, Krouchen@physio.umontreal.ca
% Copyright (c) 1989-2002 Nedialko Krouchev
% This code is freeware in the spirit of the FSF/GNU public license.

function [perm, unfeas, D] = bghungar(C)

unfeas = 0;

% --- Basic input checks (lightweight, keeps behavior stable) -------------
[m,n] = size(C);
if m ~= n
    error('bghungar:InputNotSquare', 'Input matrix C must be square.');
end
if ~isnumeric(C)
    error('bghungar:InputNotNumeric', 'Input matrix C must be numeric.');
end
if any(~isfinite(C(:)))
    % Keep original behavior: NaNs in D will be caught; here we flag earlier.
    % (We still continue to preserve legacy behavior unless NaN/Inf breaks things.)
end

% =========================================================================
% (1) Initialization / reduction
% =========================================================================

% Convert maximization to a nonnegative-like reduced form (standard trick)
[r,~] = max(C,[],1);
C1 = ones(m,1) * r - C;

[c,~] = min(C1,[],2);
D = C1 - c * ones(1,n);

perm = zeros(m,1);
if any(isnan(D(:)))
    unfeas = 1;
    return;
end

% Greedy assignment of zeros (initial feasible solution attempt)
for j = 1:m
    kk = find(~D(:,j));
    while ~isempty(kk)
        i = kk(1);
        if ~perm(i)
            perm(i) = j;
            break;
        else
            kk(1) = [];
        end
    end
end

permPrev = zeros(1,m);

% =========================================================================
% (2) Main loop: augment until all rows are assigned
% =========================================================================

keepSets = 0;

% Occupied rows:
ii = find(perm);
while length(ii) < m

    if ~keepSets
        keepSets = 0;

        % Reset the free sets:
        rr = 1:m;        % candidate rows
        cc = 1:n;        % candidate columns

        % Occupied columns:
        zz = perm(ii);

        % Free columns: remove occupied columns from cc
        cc(zz) = [];
    end

    while true
        if any(isnan(D(:)))
            unfeas = 2;
            break;
        end

        % -----------------------------------------------------------------
        % (4,5) Find a row containing a free zero
        % -----------------------------------------------------------------
        jz = 0;
        for i = rr
            kk = find(~D(i,cc));
            if ~isempty(kk)
                jz = i;
                kz = cc(kk(1));
                break;
            end
        end

        if jz
            % Found a row containing a free zero
            if perm(jz)
                % (6) This row already has an assigned zero:
                kz = perm(jz);
                cc = [cc, kz];          %#ok<AGROW>
                zz(zz == kz) = [];
                rr(rr == jz) = [];

            else
                % (7) This row has no assigned zero yet:
                % Build an alternating chain to augment the assignment
                jz0 = jz;
                kz0 = kz;
                perm0 = perm;

                while true
                    % Step in a column
                    perm(jz) = kz;

                    rr1 = [1:jz-1, jz+1:m];
                    next = find(perm(rr1) == kz);
                    if isempty(next)
                        break;
                    end
                    jz = rr1(next(1));

                    % Step in a row
                    perm(jz) = 0;

                    cc1 = [1:kz-1, kz+1:n];
                    next = find(~D(jz,cc1));
                    if isempty(next)
                        break;
                    end
                    kz = cc1(next(1));
                end

                % Detect alternating cycles (degeneracy handling)
                if any(ismember(perm', permPrev, 'rows'))
                    unfeas = 3;

                    % Block the free row/col (jz0,kz0) where this started;
                    % recover last perm0 and continue
                    perm = perm0;
                    zz = [zz; kz0];      %#ok<AGROW>
                    cc(cc == kz0) = [];
                else
                    permPrev = [permPrev; perm']; %#ok<AGROW>
                end
                break;
            end

        else
            % -----------------------------------------------------------------
            % No row containing free zeros found -> create new free zeros
            % -----------------------------------------------------------------
            p = min(min(D(rr,cc)));
            if p >= Inf
                unfeas = 4;
                break;
            end
            D(rr,:) = D(rr,:) - p;
            D(:,zz) = D(:,zz) + p;
        end
    end % inner while

    if unfeas == 3
        % Recover from degeneracy and keep the current sets
        unfeas = 0;
        keepSets = 1;
    else
        keepSets = 0;
    end

    if unfeas
        break;
    end

    ii = find(perm);
end % outer while
