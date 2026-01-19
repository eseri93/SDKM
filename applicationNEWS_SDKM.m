%% ============================================================
%  SDKM on 20 Newsgroups (sample) - Application script
%  - Loads the term-document matrix from CSV
%  - (Optional) Computes pseudo-F over a grid of (K,Q)
%  - Runs SDKM for (K,Q) = (3,2)
%  - Saves outputs to .mat
%% ============================================================

%% ----------------------------
% 1) Load data
% ----------------------------
T = readtable("data_news_sample.csv");

% First column contains terms (or IDs); the rest is the numeric matrix
terms = T(:,1);
X = T(:,2:end);
X = table2array(X);

% Set seed (legacy RNG for backward compatibility)
seed = 123;
rand('state', seed);

%% ----------------------------
% 2) Pseudo-F grid search (K,Q = 2..6)
% ----------------------------
K_list = 2:6;
Q_list = 2:6;
Rndst_grid = 10;     % smaller multistart for the grid (speed)
Stand_opt  = 'off';  % keep consistent with the final run

fprintf('\nPseudo-F grid search (20 Newsgroups sample)\n');
fprintf('Stand = %s | Rndst = %d\n\n', Stand_opt, Rndst_grid);

pF_table = NaN(length(K_list), length(Q_list));

for a = 1:length(K_list)
    for b = 1:length(Q_list)
        K = K_list(a);
        Q = Q_list(b);

        % Run SDKM for this (K,Q)
        [Vtmp, Utmp, ~, ftmp, ~] = SDKM(X, K, Q, 'Stand', Stand_opt, 'Rndst', Rndst_grid);

        % Compute pseudo-F (K and Q inferred from Utmp and Vtmp)
        pF = pseudo_F_index(X, Utmp, Vtmp);

        pF_table(a,b) = pF;

        fprintf('K=%d, Q=%d | pseudo-F=%.6g | fit=%.6g\n', K, Q, pF, ftmp);
    end
end

%% ----------------------------
% 3) Final run (K=3, Q=2)
% ----------------------------
K_final = 3;
Q_final = 2;
Rndst_final = 20;

[Vsdkm_news, Usdkm_news, Ymsdkm_news, fsdkm_news, ~] = ...
    SDKM(X, K_final, Q_final, 'Stand', Stand_opt, 'Rndst', Rndst_final);

%% ----------------------------
% 4) Quick checks / summaries
% ----------------------------
fprintf('\nFinal SDKM run (K=%d, Q=%d)\n', K_final, Q_final);
disp('Cluster sizes (units):');
disp(sum(Usdkm_news, 1));

disp('Cluster sizes (variables):');
disp(sum(Vsdkm_news, 1));

%% ----------------------------
% 5) Save outputs
% ----------------------------
save("Vsdkm_news.mat",  'Vsdkm_news');
save("Usdkm_news.mat",  'Usdkm_news');
save("Ymsdkm_news.mat", 'Ymsdkm_news');
save("fsdkm_news.mat",  'fsdkm_news');

% Optionally save pseudo-F table too
save("pF_grid_news.mat", "pF_table", "K_list", "Q_list");
