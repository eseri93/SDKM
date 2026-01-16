%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Double k-Means %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% algorithm for simultaneous partitioning of objects and variables
% tesi
function [Vdskm,Udskm,Ymdskm, fdskm,indskm]=SDKM_tesi(X, K, Q, varargin)
%

% n = number of objects
% J = number of variables
% K = number of ObjectClasses
% m = number of VariableClasses
% X = n x J (objects x variables) matrix
% U = n x K (objects x ObjectClasses) matrix
% V = J x Q (variables x variableClasses) matrix
% Rndstart = number of multistart

% centroids Xm= U'*X*V 
%
% problem: 
% max tr(X'*U*Ym*V')/sqrt(tr(X'*X)*tr(V*Ym'*U'*U*Ym*V'))
% subject to
% V binary and row stochastic
% U binary and row stochastic

%equivalent to

% max||U*Ym*V'||^2
% subject to
% V binary and row stochastic
% U binary and row stochastic

% initialization

[n,J]=size(X);
un=ones(n,1);
uk=ones(K,1);
um=ones(Q,1);
%rand('state',29);
%randn('state',29);

% 'Stats'    ->    Default value: 'on', print the statistics of the fit of the
%                  model.
%                  If 'off' Statistics are not printed (used for simulation
%                  studies)
% 'Stand'    ->    Default value 'on', standardize variables 
%                  If 'Dst' Double Standardization               
%                  If 'off' does not standardize variables 
% 
% 'Rndst'    ->    an integer values indicating the intital random starts.
%                  Default '20' thus, repeat the anaysis 20 times and retain the
%                  best solution.
% 'MaxIter'  ->    an integer value indicationg the maximum number of
%                  iterations of the algorithm
%                  Default '100'.
% 'ConvToll' ->    an arbitrary samll values indicating the convergence
%                  tollerance of the algorithm, Default '1e-9'.


%

if nargin < 2
   error('Too few inputs');
end

if ~isempty(X)
    if ~isnumeric(X)
        error('Invalid data matrix');
    end  
    if min(size(X)) == 1
    error(message('Disjoint Factor Analysis:NotEnoughData'));
end

else
    error('Empty input data matrix');
end

if ~isempty(Q)
    if isnumeric(K)
        if Q > J 
              error('The number of latent factors larger that the number of variables');
           end    
    elseif Q < 1 
              error('Invalid number of latent factors');
    end
else
    error('Empty input number of latent factors');
end

% Optional parameters   
pnames = {'Stats' 'Stand' 'Rndst' 'MaxIter' 'ConvToll'};
dflts =  { 'on'    'on'     20       100       1e-9 };
[Stats,Stand,Rndst,MaxIter,ConvToll] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%if ~isempty(eid)
%    error(sprintf('Disjoint Factor Analysis; %s',eid), emsg);
%end


% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Stats)
    if ischar(Stats)
       StatsNames = {'off', 'on'};
       js = strcmpi(Stats,StatsNames);
           if sum(js) == 0
              error(['Invalid value for the ''Statistics'' parameter: '...
                     'choices are ''on'' or ''off''.']);
           end
       Stats = StatsNames{js}; 
    else  
        error(['Invalid value for the ''Statistics'' parameter: '...
               'choices are ''on'' or ''off''.']);
    end
else 
    error(['Invalid value for the ''Statistics'' parameter: '...
           'choices are ''on'' or ''off''.']);
end
% end statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standardization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Stand)
    if ischar(Stand)
       StandNames = {'off', 'on', 'Mahalanobis','DSt'};
       js = strcmpi(Stand,StandNames);
           if sum(js) == 0
              error(['Invalid value for the ''Standardization'' parameter: '...
                     'choices are ''on'' or ''off'' or ''Mahalanobis''.']);
           end
       Stand = StandNames{js}; 
       switch Stand
           
       case 'off'
           Xs = X - ones(n,1)*mean(X);            
       case 'on'
           Xs = zscore(X,1);          
       case 'Mahalanobis'
           Jc=eye(n)-(1/n)*ones(n);
           Sx=cov(X,1);
           Xs = Jc*X*Sx^-0.5; 
       case 'DSt'
           [Xs,itt,f2]=DStand(X);
           Xs=Xs*1./sqrt(n);
       end
    else  
        error(['Invalid value for the ''standardization'' parameter: '...
               'choices are ''on'' or ''off'' or ''.']);
    end
else 
    error(['Invalid value for the ''standardization'' parameter: '...
            'choices are ''on'' or ''off'' or ''.']);
end
% end Standardization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rndst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Rndst)  
    if isnumeric(Rndst)
       if (Rndst < 0) || (Rndst > 1000) 
       error('Rndst must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Random Starts');
    end
else
    error('Invalid Number of Random Starts')
end
% end Rndst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAxIter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(MaxIter)  
    if isnumeric(MaxIter)
       if (MaxIter < 0) || (MaxIter > 1000) 
       error('MaxIter must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Max Iterations');
    end
else
    error('Invalid Number of Max Iterations')
end
% end MaxIter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ConvToll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(ConvToll)  
    if isnumeric(ConvToll)
       if (ConvToll < 0) || (ConvToll > 0.1) 
       error('ConvToll must be a value in the interval [0,0.1]');
       end
    else
       error('Invalid Convergence Tollerance');
    end
else
    error('Invalid Convergence Tollerance')
end
% end ConvToll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total deviance  
st=sum(sum(Xs.^2));


% Start the algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loop=1:Rndst
    
    V=randPU(J,Q); % generate a random partition for variables (col.s)
    U=randPU(n,K); % generate a random partition for units     (rows) 
    %su=sum(U);
    %sv=sum(V);
    %
    UC=eye(K);
    VC=eye(Q);
    KK=[1:K];
    QQ=[1:Q];
    flg=0; % increases if there is a potential null cluster
    %idu=kmeans(Xs,K,'Replicates',5);
    %U=UC(idu,:);
    %idv=kmeans(Xs',Q,'Replicates',5);
    %V=UC(idv,:);

    %[~,U,~,~]=kmeansVICHI(Xs,K,1);
    %[~,V,~,~]=kmeansVICHI((Xs'),Q,1);
    su=sum(U);
    sv=sum(V);
 
    % compute initial centroid
    dD=sqrt(diag((U'*Xs)*(Xs'*U)));
    dF=sqrt(diag((V'*Xs')*(Xs*V)));
    %dU=sqrt(diag((U'*U)));
    %dV=sqrt(diag((V'*V)));
    dU=diag((U'*U));
    dV=diag((V'*V));

    Dm1=diag(1./dD);
    Fm1=diag(1./dF);
    Um1=diag(1./dU);
    Vm1=diag(1./dV);
    Ym=normr(Um1*U'*Xs*V*Vm1);
    Xt=U*Ym*V';        
    f0=trace(Xs'*Xt)./sqrt(trace(Xs'*Xs)*trace(Xt'*Xt));

    %it=0;
    % iteration phase
            
        for it = 1:MaxIter
            % given Ym and V update U
            %Ymv=Ym*V';
            Ymv=Dm1*U'*Xs;
            %Ymv=Ym*diag(dF)*V';

            for i=1:n
                posmax=KK(U(i,:)==1); 
                maxdif=Xs(i,:)*Ymv(posmax,:)'./sqrt(Xs(i,:)*Xs(i,:)'*Ymv(posmax,:)*Ymv(posmax,:)');
                if sum(U(:,posmax))>1  
                    for k=find(KK~=posmax)
                        U(i,:)=UC(k,:); 
                        %B=U*Ym*V';
                         dif=Xs(i,:)*Ymv(k,:)'./sqrt(Xs(i,:)*Xs(i,:)'*Ymv(k,:)*Ymv(k,:)');
                        if dif > maxdif
                            maxdif=dif;
                            posmax=k;
                        end   
                    end
                    U(i,:)=UC(posmax,:);
                end
            end
             %new:
            U=balanceCluster(U);
             
            % given U and V update Ym
            dD=sqrt(diag((U'*Xs)*(Xs'*U)));
            dF=sqrt(diag((V'*Xs')*(Xs*V)));
            %dU=sqrt(diag((U'*U)));
            %dV=sqrt(diag((V'*V)));
            dU=diag((U'*U));
            dV=diag((V'*V));

            Dm1=diag(1./dD);
            Fm1=diag(1./dF);
            Um1=diag(1./dU);
            Vm1=diag(1./dV);
            Ym=normr(Um1*U'*Xs*V*Vm1);
            Xt=U*Ym*V';        
            f=trace(Xs'*Xt)./sqrt(trace(Xs'*Xs)*trace(Xt'*Xt));

            % given Ym and U update V
    
            %Ymu=U*Ym;
            Ymu=Xs*V*Fm1;
            %Ymu=U*diag(dD)*Ym;
            
            for j=1:J
                posmax=QQ(V(j,:)==1);
                maxdif=Xs(:,j)'*Ymu(:,posmax)./sqrt(Xs(:,j)'*Xs(:,j)*Ymu(:,posmax)'*Ymu(:,posmax));
                if sum(V(:,posmax))>1  
                    for i=find(QQ~=posmax)%1:Q
                        V(j,:)=VC(i,:);
                        dif=Xs(:,j)'*Ymu(:,i)./sqrt(Xs(:,j)'*Xs(:,j)*Ymu(:,i)'*Ymu(:,i));
                        if dif > maxdif
                            maxdif=dif;
                            posmax=i;
                        end   
                    end
                    V(j,:)=VC(posmax,:);
                end
            end

             %new:
            V=balanceCluster(V);

            % given U and V update Ym
            dD=sqrt(diag((U'*Xs)*(Xs'*U)));
            dF=sqrt(diag((V'*Xs')*(Xs*V)));
            %dU=sqrt(diag((U'*U)));
            %dV=sqrt(diag((V'*V)));
            dU=diag((U'*U));
            dV=diag((V'*V));

            Dm1=diag(1./dD);
            Fm1=diag(1./dF);
            Um1=diag(1./dU);
            Vm1=diag(1./dV);
            Ym=normr(Um1*U'*Xs*V*Vm1);
            Xt=U*Ym*V';        
            f=trace(Xs'*Xt)./sqrt(trace(Xs'*Xs)*trace(Xt'*Xt));
            %   
            fdif = abs(f-f0);
        
            if fdif > ConvToll 
                f0=f; 
            else
                break
                % if fdif >= 0
                %     break
                % else  
                %     [sum(U) sum(V)]
                %     %icac=[];
                %     %cac=[sum(U) sum(V)];
                %     %icac1=find(cac==1);
                %     %if isempty(icac1)  
                %         break
                %     %end
                % end
            end
        end
        disp(sprintf('DKM: Loop=%g, Explained Variance =%g, iter=%g, fdif=%g',loop,f, it,fdif))   
        if loop==1
            Vdskm=V;
            Udskm=U;
            Ymdskm=Ym;
            fdskm=f;
            loopdskm=1;
            indskm=it;
            fdifo=fdif;
        end
        if f > fdskm
            Vdskm=V;
            Udskm=U;
            Ymdskm=Ym;
            fdskm=f;
            loopdskm=loop;
            indskm=it;
            fdifo=fdif;
        end
end
% sort clusters of variables per descending order of cardinality
[~,ic]=sort(diag(Vdskm'*Vdskm), 'descend');
Vdskm=Vdskm(:,ic);
% sort clusters of objects in descending order of cardinality
[~,ic]=sort(diag(Udskm'*Udskm), 'descend');
Udskm=Udskm(:,ic);
disp(sprintf('DKM (Final): Explained Variance =%g, loopdpca=%g, iter=%g, fdif=%g',fdskm, loopdskm, indskm,fdifo)) 


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means algorithm  VICHI  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loopOtt,UOtt,fOtt,iterOtt]=kmeansVICHI(X,K,Rndstart)
%
% n = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations
%

maxiter=100;
n = size(X,1);
J = size(X,2);
epsilon=0.000001;


% initial partition U0 is given

% best in a fixed number of partitions
%seed=200;             %si può rimettere per ritrovare le soluzioni ottenute
%rand('state',seed)    %
for loop=1:Rndstart
   U0=randPU(n,K);
%
   su=sum(U0);
   %
   % given U compute Xmean (compute centroids)
   Xmean0 = diag(1./su)*U0'*X;

   for iter=1:maxiter
   %
   % given Xmean0 assign each units to the closest cluster
   %
        U=zeros(n,K);
        for i=1:n
            mindif=sum((X(i,:)-Xmean0(1,:)).^2);
            posmin=1;
            for j=2:K
                dif=sum((X(i,:)-Xmean0(j,:)).^2);
                if dif < mindif
                    mindif=dif;
                    posmin=j;
                end 
            end
            U(i,posmin)=1;
        end
   % given a partition of units 
   % i.e, given U compute Xmean (compute centroids)
   %
        su=sum(U);
        while sum(su==0)>0,
            [m,p1]=min(su);
            [m,p2]=max(su);
            ind=find(U(:,p2));
            ind=ind(1:floor(su(p2)/2));
            U(ind,p1)=1;
            U(ind,p2)=0;
            su=sum(U);
        end 
   %
   % given U compute Xmean (compute centroids)
        Xmean = diag(1./su)*U'*X;
   
   %
   % 
   % compute objective function
   %
        BB=U*Xmean-X;
        f=trace(BB'*BB);
%   
%  stopping rule
%   
        dif=sum(sum((Xmean-Xmean0).^2));
        if dif > epsilon  
            Xmean0=Xmean;
        else
            break
        end
   end
   if loop==1
        UOtt=U;
        fOtt=f;
        loopOtt=1;
        iterOtt=1;
        XmeanOtt=Xmean;
   end
   if f < fOtt
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double standardization          %
% July 2016                       %
% M Vichi                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xs,it,f]=DStand(X)

% Example %%%%%%%
%X= [0.1182 0.7069 0.4145;
%    0.9884 0.9995 0.4648;
%    0.5400 0.2878 0.7640];


eps=1e-14;
[n,J]=size(X);

un=ones(n,1);
uJ=ones(J,1);
Xso=X;
for it=1:30000
    % Standardize Columns 
    % column centring
    X1=zscore(Xso,1);
    % Standardize Rows
    Xs=(zscore(X1',1))';
    f=norm(Xs-Xso);
    if f < eps
        break
    end
    Xso=Xs;
end
[it f];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Balance Clusters          %
% ????                    %
% M Vichi                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A]=balanceCluster(A)
    sv=sum(A);
    while any(sv==0)
        [~,p1]=min(sv);
        [~,p2]=max(sv);
        ind = find(A(:,p2));              
        ind=ind(1:floor(sv(p2)/2));
        A(ind,p1)=1;
        A(ind,p2)=0;
        sv=sum(A);
    end

