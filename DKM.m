%%%%%%%%%%%%%%%%%%
% Double k-Means %
%%%%%%%%%%%%%%%%%%
%
% algorithm for simultaneous partitioning of objects and variables
% Maurizio Vichi October 2012
%
function [Vdkm,Udkm,Ymdkm, fdkm,indkm]=DKM(X, K, Q, Rndstart)
%
% n = number of objects
% J = number of variables
% K = number of ObjectClasses
% m = number of VariableClasses
% X = n x J (objects x variables) matrix
% U = n x K (objects x ObjectClasses) matrix
% V = J x Q (variables x variableClasses) matrix
%
% model X = U*Ym*V' + E
%
% problem: min||X-U*Ym*V'||^2

%equivalent to

% max||U*Ym*V'||^2
% subject to
% V binary and row stochastic
% U binary and row stochastic

% initialization
%
maxiter=100;
% convergence tollerance
eps=0.0000000001;


[n,J]=size(X);

% centring matrix
Jm=eye(n)-(1/n)*ones(n);

% compute var-covar matrix 
S=(1/n)*X'*Jm*X;

% Standardize data
Xs=Jm*X*diag(diag(S))^-0.5;


st=sum(sum(Xs.^2));

un=ones(n,1);
uk=ones(K,1);
um=ones(Q,1);


for loop=1:Rndstart
    V=randPU(J,Q);
    U=randPU(n,K);
    su=sum(U);
    sv=sum(V);
    Ym = diag(1./su)*U'*Xs*V*diag(1./sv);
    B=U*Ym*V';
    f0=trace(B'*B)/st;
    it=0;
    % iteration phase

        fdif=2*eps;
        while fdif > eps | it>=maxiter,
            it=it+1;
           
            % given Xm and V update U
            U=zeros(n,K);
            Ymv=Ym*V';
            for i=1:n
                mindif=sum((Xs(i,:)-Ymv(1,:)).^2);
                posmin=1;
                for k=2:K
                    dif=sum((Xs(i,:)-Ymv(k,:)).^2);
                    if dif < mindif
                        mindif=dif;
                        posmin=k;
                    end 
                end
                U(i,posmin)=1;
            end
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
            
            % given U and V updata Xm
      
            su=sum(U);
            Ym = diag(1./su)*U'*Xs*V*diag(1./sv);  

            % diven Xm and U updata V
    
            V=zeros(J,Q);
            Ymu=U*Ym;
            for j=1:J
                mindif=sum((Xs(:,j)-Ymu(:,1)).^2);
                posmin=1;
                for i=2:Q
                    dif=sum((Xs(:,j)-Ymu(:,i)).^2);
                    if dif < mindif
                        mindif=dif;
                        posmin=i;
                    end 
                end
                V(j,posmin)=1;
            end

            sv=sum(V);
            while sum(sv==0)>0,
                [m,p1]=min(sv);
                [m,p2]=max(sv);
                ind=find(V(:,p2));
                ind=ind(1:floor(sv(p2)/2));
                V(ind,p1)=1;
                V(ind,p2)=0;
                sv=sum(V);
            end 

            % given U and V updata Xm
      
            sv=sum(V);
            Ym = diag(1./su)*U'*Xs*V*diag(1./sv);

            B=U*Ym*V';
            f=trace(B'*B)/st;
%   
            fdif = f-f0;
        
            if fdif > eps 
                f0=f; 
            else
                break
            end
        end
        disp(sprintf('DKM: Loop=%g, Explained Variance =%g, iter=%g, fdif=%g',loop,f, it,fdif))   
        if loop==1
            Vdkm=V;
            Udkm=U;
            Ymdkm=Ym;
            fdkm=f;
            loopdkm=1;
            indkm=it;
            fdifo=fdif;
        end
        if f > fdkm
            Vdkm=V;
            Udkm=U;
            Ymdkm=Ym;
            fdkm=f;
            loopdkm=loop;
            indkm=it;
            fdifo=fdif;
        end
end
% sort clusters of variables per descending order of cardinality
[c,ic]=sort(diag(Vdkm'*Vdkm), 'descend');
Vdkm=Vdkm(:,ic);
% sort clusters of objects in descending order of cardinality
[c,ic]=sort(diag(Udkm'*Udkm), 'descend');
Udkm=Udkm(:,ic);

disp(sprintf('DKM (Final): Explained Variance =%g, loopdpca=%g, iter=%g, fdif=%g',fdkm, loopdkm, indkm,fdifo)) 