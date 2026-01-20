%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose_RndStarts.m
%
% Empirical choice of the number of random starts (Rndst) in SDKM to reduce
% the risk of local maxima. The script:
%   1) Generates one synthetic dataset from GenDataSDKM
%   2) Computes the "true" objective value (using the true partitions)
%   3) Runs SDKM repeatedly for different Rndst values
%   4) Counts how often the obtained objective is lower than the true one
%
% NOTE:
% - This is a diagnostic script for tuning Rndst (not required for running
%   applications).
% - Requires: GenDataSDKM.m, SDKM.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ----------------------------
% 0) Reproducibility
% ----------------------------
seed = 123;
rng(seed, "twister");        % modern RNG (recommended)
% If you prefer legacy compatibility, you could use:
% rand('state', seed);

%% ----------------------------
% 1) Synthetic data generation
% ----------------------------
n      = 100;     % number of units (objects)
J      = 30;      % number of variables
K_true = 3;       % true number of row clusters
Q_true = 2;       % true number of column clusters
er1    = 1.5;     % centroid heterogeneity
er2    = 1.5;     % noise / within-cluster heterogeneity

% Generate one dataset with known ground-truth partitions
[X, Ym_true, U_true, V_true] = GenDataSDKM(n, K_true, J, Q_true, er1, er2);

%% ----------------------------
% 2) Compute the "true" objective value
% ----------------------------
% Standardize X as done in SDKM when 'Stand' is 'on'
Xs = zscore(X, 1);

% Reconstruction using the TRUE co-clustering structure
Xt = U_true * Ym_true * V_true';

% Cosine-like objective (same structure used in SDKM output fsdkm)
% f = <Xs, Xt> / (||Xs||_F * ||Xt||_F)
num   = trace(Xs' * Xt);
denom = sqrt(trace(Xs' * Xs) * trace(Xt' * Xt));
f_true = num / denom;

%% ----------------------------
% 3) Run SDKM for different numbers of random starts
% ----------------------------
rndS_values = [1 5 10 20 30 40 50 70 100];  % candidate Rndst values
num_runs    = 100;                          % repeated runs per setting

% Store SDKM objective values: rows = rndS_values, cols = runs
obj_f = zeros(length(rndS_values), num_runs);

for r = 1:num_runs
    for i = 1:length(rndS_values)
        rndS = rndS_values(i);

        % Run SDKM
        % (K_true, Q_true fixed here because we focus on local maxima vs starts)
        [~, ~, ~, fsdkm, ~] = SDKM(X, K_true, Q_true, 'Stand', 'on', 'Rndst', rndS);

        obj_f(i, r) = fsdkm;
    end
end

%% ----------------------------
% 4) Summarize: how often SDKM is below the "true" objective
% ----------------------------
below_true_counts = sum(obj_f < f_true, 2);          % count per rndS setting
below_true_table  = table(rndS_values(:), below_true_counts, ...
    'VariableNames', {'Rndst', 'NumRunsBelowTrue'});

disp("Counts of runs with fsdkm < f_true (lower suggests local maxima / suboptimal):");
disp(below_true_table);

%% ----------------------------
% 5) Save output
% ----------------------------
outdir = "Output";
if ~exist(outdir, "dir"), mkdir(outdir); end

save(fullfile(outdir, "scelta_RndStarts.mat"), ...
     "seed","n","J","K_true","Q_true","er1","er2", ...
     "rndS_values","num_runs","obj_f","f_true","below_true_counts");

fprintf("Saved results to %s\n", fullfile(outdir, "scelta_RndStarts.mat"));

