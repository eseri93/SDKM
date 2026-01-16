%% different scenarios simulation study%%
%% high level
%% data generation
% Parameters for data generation
seed=123;
rand('state',seed);    %
n = 100; % Number of units (objects)
J = 30; % Number of variables
K_true = 3; % True number of clusters for units
Q_true = 2; % True number of clusters for variables
er1 = 0.1; % Heterogeneity of centroids
er2 = 2; % Heterogeneity of clusters
%er1=0.1 and er2=2 permette di aver mean ari of U 0.59 and mean ari of
%V=0.71
% er1=0.1, er2=1.5 permette di avere % di local maxima=0 con Random Start
% 20.
% Generate synthetic data with known clusters for each run
[X, Ym_true, UdkmO, VdkmO] = GenDataSDKM(n, K_true, J, Q_true, er1, er2);
%% compute ari for comparing partitions and rmse for comparing centroids matrix
num_runs = 100; % Number of runs
ari_u= zeros(num_runs,1);
ari_v=zeros(num_runs,1);
centroids = cell(num_runs, 1);
rmse_v=ones(num_runs,1);
for j=1:num_runs
[Vsdkm, Usdkm, Ymsdkm, fsdkm, ~] = SDKM_tesi(X, K_true, Q_true, 'Stand', 'on', 'Rndst', 20);
ari_u(j)=mrand(Usdkm' * UdkmO); % ARI for unit clusters
ari_v(j)=mrand(Vsdkm' * VdkmO); % ARI for variables clusters
centroids{j}=Ymsdkm;
rmse_v(j)=rmse(Ymsdkm(:),Ym_true(:));
end
mean(ari_u)
mean(ari_v)
median(ari_u)
median(ari_v)
%la funzione rmse fa la seguente operazione
% (y - yhat) % Errors
% (y - yhat).^2 % Squared Error
% mean((y - yhat).^2) % Mean Squared Error
% RMSE = sqrt(mean((y - yhat).^2)); % Root Mean Squared Error
% e restiruisce l'RMSE di colonna
% per avere l'RMSE totale metto (:) dopo le matrici (in questo modo le
% rendo vettori). E ottendo un unico valore di RMSE
mean(rmse_v)
median(rmse_v)


%% low level
%% data generation
% Parameters for data generation
seed=123;
rand('state',seed);    %
n = 100; % Number of units (objects)
J = 30; % Number of variables
K_true = 3; % True number of clusters for units
Q_true = 2; % True number of clusters for variables
er1 = 0.0001; % Heterogeneity of centroids
er2 = 0.001;
%er1=0.1 and er2=2 permette di aver mean ari of U 0.59 and mean ari of
%V=0.71
% er1=0.1, er2=1.5 permette di avere % di local maxima=0 con Random Start
% 20.
% Generate synthetic data with known clusters for each run
[X, Ym_true, UdkmO, VdkmO] = GenDataSDKM(n, K_true, J, Q_true, er1, er2);
%% compute ari for comparing partitions and rmse for comparing centroids matrix
num_runs = 100; % Number of runs
ari_u= zeros(num_runs,1);
ari_v=zeros(num_runs,1);
centroids = cell(num_runs, 1);
rmse_v=ones(num_runs,1);
for j=1:num_runs
[Vsdkm, Usdkm, Ymsdkm, fsdkm, ~] = SDKM_tesi(X, K_true, Q_true, 'Stand', 'on', 'Rndst', 20);
ari_u(j)=mrand(Usdkm' * UdkmO); % ARI for unit clusters
ari_v(j)=mrand(Vsdkm' * VdkmO); % ARI for variables clusters
centroids{j}=Ymsdkm;
rmse_v(j)=rmse(Ymsdkm(:),Ym_true(:));
end
mean(ari_u)
mean(ari_v)
median(ari_u)
median(ari_v)
%la funzione rmse fa la seguente operazione
% (y - yhat) % Errors
% (y - yhat).^2 % Squared Error
% mean((y - yhat).^2) % Mean Squared Error
% RMSE = sqrt(mean((y - yhat).^2)); % Root Mean Squared Error
% e restiruisce l'RMSE di colonna
% per avere l'RMSE totale metto (:) dopo le matrici (in questo modo le
% rendo vettori). E ottendo un unico valore di RMSE
mean(rmse_v)
median(rmse_v)
