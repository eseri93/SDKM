%% ============================================================
%  Application: SDKM on U.S. Presidential Inaugural Addresses
%  - Load data from CSV
%  - Compute pseudo-F for K,Q in {2,...,6}
%  - Run SDKM for K=3, Q=2 and save results
% ============================================================

%% Read data
Xtbl = readtable("data_text.csv");

% Store terms (first column) and extract numeric matrix (remaining columns)
terms = Xtbl(:,1); %#ok<NASGU>  % kept for possible downstream use
X = Xtbl(:,2:end);
X = table2array(X);

% Reproducibility (legacy + modern RNG)
seed = 123;
rand('state', seed);     %#ok<RAND>
rng(seed, 'twister');

[n, J] = size(X);

%% ------------------------------------------------------------
%  Pseudo-F grid search over K and Q in {2,...,6}
%  NOTE: SDKM with 'Stand'='off' centers columns internally.
%  To keep pseudo-F consistent with that setting, we compute it
%  on the column-centered matrix here.
%% ------------------------------------------------------------
Xc = bsxfun(@minus, X, mean(X, 1));  % column-centered version of X

K_list = 2:6;
Q_list = 2:6;

pF_grid = NaN(numel(K_list), numel(Q_list));
f_grid  = NaN(numel(K_list), numel(Q_list));  % optional: store explained variance too

% You can increase Rndst_grid if you want a more thorough search
Rndst_grid = 10;

fprintf('\n============================================================\n');
fprintf('Pseudo-F grid search (K,Q in 2..6) using SDKM solutions\n');
fprintf('Data size: n=%d, J=%d | Rndst_grid=%d | Stand=off\n', n, J, Rndst_grid);
fprintf('============================================================\n\n');

for a = 1:numel(K_list)
    K = K_list(a);
    for b = 1:numel(Q_list)
        Q = Q_list(b);

        % Skip infeasible configurations
        if (K > n) || (Q > J)
            continue;
        end

        % Run SDKM for this (K,Q)
        [Vtmp, Utmp, ~, ftmp, ~] = SDKM(X, K, Q, 'Stand', 'off', 'Rndst', Rndst_grid);

        % Compute pseudo-F for the obtained (U,V) on centered X (consistent with Stand='off')
        pF = pseudo_F_index(Xc, Utmp, Vtmp);

        pF_grid(a,b) = pF;
        f_grid(a,b)  = ftmp;

        fprintf('K=%d, Q=%d | pseudo-F = %.6g | explained var = %.6g\n', K, Q, pF, ftmp);
    end
end

% Print pseudo-F table
rowNames = compose('K_%d', K_list);
colNames = compose('Q_%d', Q_list);
pF_table = array2table(pF_grid, 'RowNames', rowNames, 'VariableNames', colNames);

fprintf('\nPseudo-F grid (rows=K, cols=Q):\n');
disp(pF_table);

% Report best (K,Q) by pseudo-F (global maximum over computed grid)
[bestVal, idxMax] = max(pF_grid(:));
[bestA, bestB] = ind2sub(size(pF_grid), idxMax);
bestK = K_list(bestA);
bestQ = Q_list(bestB);

fprintf('Best (global) pseudo-F in the grid: K=%d, Q=%d with pseudo-F=%.6g\n\n', bestK, bestQ, bestVal);

%% ------------------------------------------------------------
%  Run SDKM with K=3 and Q=2 (as in the original script)
%% ------------------------------------------------------------
fprintf('============================================================\n');
fprintf('Running final SDKM with K=3, Q=2 (Rndst=20, Stand=off)\n');
fprintf('============================================================\n\n');

[Vsdkm, Usdkm, Ymsdkm, fsdkm, ~] = SDKM(X, 3, 2, 'Stand', 'off', 'Rndst', 20);

%% Results analysis
sum(Usdkm)
sum(Vsdkm)

%% Save results
save("Vsdkm.mat",  'Vsdkm')
save("Usdkm.mat",  'Usdkm')
save("Ymsdkm.mat", 'Ymsdkm')
save("fsdkm.mat",  'fsdkm')

% Optional: save pseudo-F grid too
save("pseudoF_grid.mat", "pF_grid", "K_list", "Q_list", "f_grid", "Rndst_grid")
