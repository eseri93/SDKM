%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Generator Spherical Double K-means %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = number of objects
% K = number of clusters of objects
% J = number of variables 
% Q = number of cluster of variables (1<Q<min{K-1,J})
% er1: heterogeneity of centroids
% er2: heterogeneity of clusters

function [X, Ym,UdkmO,VdkmO]=GenDataSDKM(n,K,J,Q,er1,er2);
%rand('state',34);
%randn('state',34);

% generate random partitions:
UdkmO=randPU(n,K);
VdkmO=randPU(J,Q);

% sort clusters by membership of V
vV=VdkmO*[1:Q]';
[~,ic]=sort(vV);
VdkmO=VdkmO(ic,:);
% sort clusters by membership of U;
vU=UdkmO*[1:K]';
[~,ic]=sort(vU);
UdkmO=UdkmO(ic,:);

D=10*(ones(K)-eye(K));
Jc=eye(K)-(1/K)*ones(K);
cDc=-0.5*Jc*D*Jc;
[EiVe,EiVa]=eigs(cDc,Q, 'lr');
Ym=EiVe*EiVa.^0.5;

% add an error to extreme points in order to modify their distance and consequently their isolation
%
Ym=Ym+randn(K,Q)*er1;   %now almost equidistant 

%normalize by row:
Ym=normr(Ym);

% data generation
X=UdkmO*Ym*VdkmO'+randn(n,J)*er2; 


% display a Data Matrix
imagesc(X)
colormap(flipud(hot));
colorbar

function [U]=randPU(n,c)

% generates a random partition of n objects in c classes with non-empty
% classes
%
% n = number of objects
% c = number of classes
%
U=zeros(n,c);
U(1:c,:)=eye(c);

U(c+1:n,1)=1;
for i=c+1:n
    U(i,[1:c])=U(i,randperm(c));
end
U(:,:)=U(randperm(n),:);






    


