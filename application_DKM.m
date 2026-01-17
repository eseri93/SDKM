%% read data
X=readtable("data_text.csv");
%remove first column
terms=X(:,1);
X=X(:,2:end);
X=table2array(X);
seed=123;
rand('state',seed);  

%% apply DKM function
[Vdkm,Udkm,Ymdkm, fdkm,indkm]=DKM(X, 3, 2, 20);

%% results analysis
sum(Udkm)
sum(Vdkm)

%% save results
save("Vdkm.mat",'Vdkm')
save("Udkm.mat",'Udkm')
save("Ymdkm.mat",'Ymdkm')
save("fdkm.mat",'fdkm')