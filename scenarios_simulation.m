%% scenarios_simulation.m
% Simulation study for SDKM under different noise / centroid-perturbation levels.
%
% For each eps in eps_list:
%   1) (Optional) Export simulated datasets to .mat for R / vMF comparison
%   2) Run SDKM num_runs times and compute:
%        - ARI for row clusters (U) and column clusters (V)
%        - RMSE / NRMSE for centroid matrix Ym after label switching
%
% Dependencies (must be on MATLAB path):
%   - SDKM.m
%   - GenDataSDKM.m
%   - mrand.m        (Adjusted Rand Index function you use)
%   - bghungar.m     (Hungarian algorithm for label switching)

clear; clc;

%% ----------------------------
% 0) Global settings
% ----------------------------
outdir_results = fullfile('Output');         % where simulation metrics are saved
outdir_simdata = fullfile('OutputSimData');  % where exported datasets are saved

if ~exist(outdir_results, 'dir'), mkdir(outdir_results); end
if ~exist(outdir_simdata, 'dir'), mkdir(outdir_simdata); end

% Data-generation parameters
n      = 100;
J      = 30;
K_true = 3;
Q_true = 2;

% Simulation grid (centroid perturbation and noise are both set to eps)
eps_list  = [0.10 0.35 0.50 0.75 0.90 1.10 1.35 1.50 1.75 2.00];

% Monte Carlo settings
num_runs  = 500;
base_seed = 123;  % base seed (we will deterministically derive per eps/run seeds)

% SDKM options (keep consistent across eps)
standOpt = 'off';   % IMPORTANT: 'off' matches your centroid RMSE computed on raw X
rndst    = 20;

% Optional: export datasets for vMF comparison (R.matlab friendly)
do_export_for_vmf = true;

%% ----------------------------
% 1) Main loop over epsilon values
% ----------------------------
summary = struct();
summary.eps = eps_list(:);

mean_ari_u   = zeros(numel(eps_list),1);
median_ari_u = zeros(numel(eps_list),1);
mean_ari_v   = zeros(numel(eps_list),1);
median_ari_v = zeros(numel(eps_list),1);

mean_rmse    = zeros(numel(eps_list),1);
median_rmse  = zeros(numel(eps_list),1);

mean_nrmse_r = zeros(numel(eps_list),1);
median_nrmse_r = zeros(numel(eps_list),1);

mean_nrmse_n = zeros(numel(eps_list),1);
median_nrmse_n = zeros(numel(eps_list),1);

for eIdx = 1:numel(eps_list)
    eps_val = eps_list(eIdx);
    fprintf('\n=== eps = %.2f (%d/%d) ===\n', eps_val, eIdx, numel(eps_list));

    % 1A) Optional: export data stacks for R/vMF
    if do_export_for_vmf
        export_simdata_for_vmf(outdir_simdata, eps_val, n, J, K_true, Q_true, num_runs, base_seed, eIdx);
    end

    % 1B) Run simulation for this epsilon
    res = run_sdkm_simulation(eps_val, n, J, K_true, Q_true, num_runs, base_seed, eIdx, standOpt, rndst);

    % Save results for this eps
    tag = epsTag(eps_val);
    fout = fullfile(outdir_results, sprintf('results_sim_eps_%s.mat', tag));
    save(fout, '-struct', 'res');
    fprintf('Saved simulation results: %s\n', fout);

    % Store summaries for a quick table
    mean_ari_u(eIdx)   = mean(res.ari_u);
    median_ari_u(eIdx) = median(res.ari_u);
    mean_ari_v(eIdx)   = mean(res.ari_v);
    median_ari_v(eIdx) = median(res.ari_v);

    mean_rmse(eIdx)    = mean(res.rmse_Y);
    median_rmse(eIdx)  = median(res.rmse_Y);

    mean_nrmse_r(eIdx)   = mean(res.nrmse_range);
    median_nrmse_r(eIdx) = median(res.nrmse_range);

    mean_nrmse_n(eIdx)   = mean(res.nrmse_norm);
    median_nrmse_n(eIdx) = median(res.nrmse_norm);
end

%% ----------------------------
% 2) Print summary table (console)
% ----------------------------
fprintf('\n\n=== SUMMARY (num_runs = %d, Stand = %s, Rndst = %d) ===\n', num_runs, standOpt, rndst);
T = table( ...
    summary.eps, ...
    mean_ari_u, median_ari_u, ...
    mean_ari_v, median_ari_v, ...
    mean_rmse,  median_rmse, ...
    mean_nrmse_r, median_nrmse_r, ...
    mean_nrmse_n, median_nrmse_n, ...
    'VariableNames', {'eps','meanARI_U','medianARI_U','meanARI_V','medianARI_V','meanRMSE_Y','medianRMSE_Y','meanNRMSE_range','medianNRMSE_range','meanNRMSE_norm','medianNRMSE_norm'} ...
);
disp(T);

% Save summary table too
save(fullfile(outdir_results, 'summary_table.mat'), 'T');


%% ========================================================================
% Local functions (must be at the end of a script)
% ========================================================================

function res = run_sdkm_simulation(eps_val, n, J, K_true, Q_true, num_runs, base_seed, eIdx, standOpt, rndst)
%RUN_SDKM_SIMULATION  Run SDKM on num_runs datasets and compute metrics.

ari_u       = zeros(num_runs,1);
ari_v       = zeros(num_runs,1);
rmse_Y      = zeros(num_runs,1);
nrmse_range = zeros(num_runs,1);
nrmse_norm  = zeros(num_runs,1);

for r = 1:num_runs
    % Deterministic per-(eps,run) seed
    rng(makeSeed(base_seed, eIdx, r), 'twister');

    % Generate synthetic data
    [X, Ym_true, U_true, V_true] = GenDataSDKM(n, K_true, J, Q_true, eps_val, eps_val);

    % Run SDKM
    [V_hat, U_hat, ~, f_hat, it_hat] = SDKM(X, K_true, Q_true, 'Stand', standOpt, 'Rndst', rndst);

    % ARI for partitions (as in your original code)
    ari_u(r) = mrand(U_hat' * U_true);
    ari_v(r) = mrand(V_hat' * V_true);

    % Label switching correction (Hungarian on membership mismatch)
    U_ok = fixLabelSwitching(U_hat, U_true);
    V_ok = fixLabelSwitching(V_hat, V_true);

    % Centroid estimate in the ORIGINAL model space:
    % X ≈ U * Y * V'
    Y_hat = estimateY_raw(X, U_ok, V_ok);

    % Store centroid errors
    rmse_Y(r) = rmseScalar(Y_hat(:), Ym_true(:));

    denom_range = (max(Ym_true(:)) - min(Ym_true(:)));
    if denom_range == 0, denom_range = 1; end
    nrmse_range(r) = rmse_Y(r) / denom_range;

    denom_norm = norm(Ym_true(:) - mean(Ym_true(:)));
    if denom_norm == 0, denom_norm = 1; end
    nrmse_norm(r) = rmse_Y(r) / denom_norm;

    % Optional progress every 50 runs
    if mod(r,50) == 0
        fprintf('  run %d/%d | f=%.4f | it=%d | ARI(U)=%.3f | ARI(V)=%.3f\n', ...
            r, num_runs, f_hat, it_hat, ari_u(r), ari_v(r));
    end
end

% Output struct saved to file
res = struct();
res.eps_val      = eps_val;
res.n            = n;
res.J            = J;
res.K_true       = K_true;
res.Q_true       = Q_true;
res.num_runs     = num_runs;
res.base_seed    = base_seed;
res.standOpt     = standOpt;
res.rndst        = rndst;

res.ari_u        = ari_u;
res.ari_v        = ari_v;
res.rmse_Y       = rmse_Y;
res.nrmse_range  = nrmse_range;
res.nrmse_norm   = nrmse_norm;
end


function export_simdata_for_vmf(outdir, eps_val, n, J, K_true, Q_true, num_runs, base_seed, eIdx)
%EXPORT_SIMDATA_FOR_VMF  Save a stack of simulated datasets for R/vMF analyses.
% Saves: X_stack (n x J x num_runs), U_true_stack (n x K x num_runs), V_true_stack (J x Q x num_runs)
% Uses -v7 for compatibility with R.matlab.

X_stack      = zeros(n, J, num_runs, 'double');
U_true_stack = zeros(n, K_true, num_runs, 'uint8');
V_true_stack = zeros(J, Q_true, num_runs, 'uint8');

for r = 1:num_runs
    rng(makeSeed(base_seed, eIdx, r), 'twister');
    [X, ~, U_true, V_true] = GenDataSDKM(n, K_true, J, Q_true, eps_val, eps_val);

    X_stack(:,:,r)      = X;
    U_true_stack(:,:,r) = uint8(U_true);
    V_true_stack(:,:,r) = uint8(V_true);
end

params = struct('n',n,'J',J,'K_true',K_true,'Q_true',Q_true, ...
                'er1',eps_val,'er2',eps_val,'num_runs',num_runs,'base_seed',base_seed,'eps_index',eIdx);

tag = epsTag(eps_val);
fname = fullfile(outdir, sprintf('simdata_eps_%s.mat', tag));

save(fname, 'X_stack', 'U_true_stack', 'V_true_stack', 'params', '-v7');
fprintf('Saved simdata stack: %s\n', fname);
end


function A_ok = fixLabelSwitching(A_hat, A_true)
%FIXLABELSWITCHING  Reorder columns of A_hat to best match A_true using Hungarian algorithm.

c = size(A_true,2);
C = zeros(c,c);
for k = 1:c
    for g = 1:c
        C(k,g) = sum((A_true(:,k) - A_hat(:,g)).^2);
    end
end

% bghungar expects a profit matrix -> we pass -C
perm = bghungar(-C);
A_ok = A_hat(:, perm);
end


function Y = estimateY_raw(X, U, V)
%ESTIMATEY_RAW  Estimate Y from X, U, V in the model X ≈ U * Y * V'.
% Uses the same algebra as your original code, but avoids toolbox normr.

dU = diag(U' * U);
dV = diag(V' * V);

% Avoid division by zero (should not happen if clusters are non-empty)
dU(dU == 0) = 1;
dV(dV == 0) = 1;

Um1 = diag(1 ./ dU);
Vm1 = diag(1 ./ dV);

Y = Um1 * (U' * X * V) * Vm1;
Y = rowNormalize(Y);
end


function X = rowNormalize(X)
%ROWNORMALIZE  Normalize each row to have unit L2 norm (safe for zero rows).
nr = sqrt(sum(X.^2, 2));
nr(nr == 0) = 1;
X = X ./ nr;
end


function r = rmseScalar(yhat, y)
%RMSESCALAR  Scalar RMSE between two vectors.
d = yhat - y;
r = sqrt(mean(d.^2));
end


function s = epsTag(eps_val)
%EPS_TAG  Convert epsilon to a filename-safe tag like 0p10, 1p35, etc.
s = strrep(sprintf('%.2f', eps_val), '.', 'p');
end


function seed = makeSeed(base_seed, eps_index, run_index)
%MAKESEED  Deterministic per-(eps,run) seed.
seed = base_seed + 100000 * eps_index + run_index;
end
