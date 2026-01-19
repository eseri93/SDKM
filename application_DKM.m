%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application script: DKM on U.S. Presidential Inaugural Data
% - Loads the term-document matrix from data_text.csv
% - (Optional) computes pseudo-F over a grid of (K,Q)
% - Runs DKM with K=3, Q=2 and saves outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------------
% 0) Load data
% ----------------------------
T = readtable("data_text.csv");

% Keep terms (first column) for possible downstream interpretation
terms = T(:,1);

% Remove first column and convert to numeric matrix (n x J)
X = table2array(T(:,2:end));

% Set RNG for reproducibility
seed = 123;
rng(seed, 'twister');

%% ----------------------------
% 1) OPTIONAL: pseudo-F grid search (K,Q = 2..6)
% ----------------------------
do_grid_search = true;   % set to false if you want to skip this block

if do_grid_search
    % Center X (recommended for pseudo-F comparability)
    Xc = X - mean(X, 1);

    K_list = 2:6;
    Q_list = 2:6;
    Rndstart_grid = 20;

    results = []; % rows: [K Q fit pF]

    fprintf('\nPseudo-F grid search for DKM (K,Q in 2..6)\n');
    fprintf('-----------------------------------------------------------\n');
    fprintf('%4s %4s %14s %14s\n', 'K', 'Q', 'fit(fdkm)', 'pF');

    for K = K_list
        for Q = Q_list
            if K > size(X,1) || Q > size(X,2)
                continue;
            end

            [Vtmp, Utmp, ~, ftmp, ~] = DKM(X, K, Q, Rndstart_grid);

            % Pseudo-F on centered data
            pF = pseudo_F_index(Xc, Utmp, Vtmp);

            results = [results; K, Q, ftmp, pF]; %#ok<AGROW>
            fprintf('%4d %4d %14.6g %14.6g\n', K, Q, ftmp, pF);
        end
    end

    % Sort by pF descending and display best few
    [~, ord] = sort(results(:,4), 'descend');
    topN = min(10, size(results,1));

    fprintf('\nTop %d solutions by pseudo-F:\n', topN);
    fprintf('%4s %4s %14s %14s\n', 'K', 'Q', 'fit(fdkm)', 'pF');
    for i = 1:topN
        r = results(ord(i),:);
        fprintf('%4d %4d %14.6g %14.6g\n', r(1), r(2), r(3), r(4));
    end
    fprintf('-----------------------------------------------------------\n\n');
end

%% ----------------------------
% 2) Run DKM (final choice)
% ----------------------------
K_final = 3;
Q_final = 2;
Rndstart_final = 20;

[Vdkm, Udkm, Ymdkm, fdkm, indkm] = DKM(X, K_final, Q_final, Rndstart_final);

%% ----------------------------
% 3) Quick checks
% ----------------------------
disp('Cluster sizes (rows/U):');
disp(sum(Udkm, 1));

disp('Cluster sizes (cols/V):');
disp(sum(Vdkm, 1));

fprintf('Final DKM fit (explained variance proxy): %.6g\n', fdkm);
fprintf('Iterations used (best run): %d\n', indkm);

%% ----------------------------
% 4) Save outputs
% ----------------------------
save("Vdkm.mat",  "Vdkm");
save("Udkm.mat",  "Udkm");
save("Ymdkm.mat", "Ymdkm");
save("fdkm.mat",  "fdkm");
save("indkm.mat", "indkm");
