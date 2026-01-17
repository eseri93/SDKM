%% read data
X=readtable("data_text.csv");
%remove first column
terms=X(:,1);
X=X(:,2:end);
X=table2array(X);
seed=123;
rand('state',seed);    
%% apply SDKM function
[Vsdkm, Usdkm, Ymsdkm, fsdkm, ~] = SDKM_tesi(X, 3, 2, 'Stand', 'off', 'Rndst', 20);
%% results analysis
sum(Usdkm)
sum(Vsdkm)
%% save results
save("Vsdkm.mat",'Vsdkm')
save("Usdkm.mat",'Usdkm')
save("Ymsdkm.mat",'Ymsdkm')
save("fsdkm.mat",'fsdkm')