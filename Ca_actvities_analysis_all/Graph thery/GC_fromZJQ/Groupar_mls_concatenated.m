function [GC_nonzero,GCstrength,A,pF,pWald]=Groupar_mls_concatenated(TS,T0,p,Netlabel) 
% Input
% 'TS'       concatenated time course of different subjects or sessions encoded in a matrix
%            of size (T=total number of time points = T0 multiply by number of subjects/sessions
%            T0 = length of time series of one subject/session,k=number of variables/ROIs)
% 'p'             order of the AR model to be identified
% Y = A*Z+E,    Eq. (1)

% where: Y = [TS(p+1,:)' TS(p+2,:)' ... TS(T,:)'] is a matrix of size (k,T-p)
%
%        A = [A_0 A_1 ... A_p] is a matrix of size (k,k*p+1) that
%        gathers unknown parameters of the model
%
%            |    1             1           ...        1     |
%            | TS(p,:)'     TS(p+1,:)'      ...    TS(T-1,:)'|
%            |TS(p-1,:)'     TS(p,:)'       ...    TS(T-2,:)'|
%        Z = |    .             .          .           .     |
%            |    .             .            .         .     |
%            |    .             .              .       .     |
%            | TS(1,:)'      TS(2,:)'       ...    TS(T-p,:)'|
%
%        is a matrix of size (k*p+1,T-p) that is directly built from the input TS.
%
%        E = [e(p+1) e(p+2) ... e(T)] is a matrix of size (k,T-p)
%        containing the residuals of the multivariate AR model.
% Output
%
% 'GC_nonzero'    Logical Matrix consisting of 0-1.Where GC_nonzero(i,j)=1 
%                 represents a causal link from j to i
% 'GCstrength'    Matrix representing the Granger causality strenth from
%                 column index to raw index
% 'A'             Matrix containing AR model parameters
% 'pF'            Matrix of p values of F test
% 'pWald'         Matrix of p values of Wald test
if exist('p')~=1 p=1;end
T   = size(TS,1);%length of time series.T must be an integral multiple of T0
k   = size(TS,2);%number of ROIs
n   = T/T0;%number of subjects
% Normalize the  time series
for i=1:n TS(T0*(i-1)+1:T0*i,:)=zscore(TS(T0*(i-1)+1:T0*i,:));end
Y=[];
for i=1:n
    Y=[Y,TS(T0*(i-1)+p+1:T0*i,:)'];
end

%concatenate the predictors
Z(1,:)=ones(1,n*(T0-p));
for j=1:p
    for t=1:n
    Z(2+k*(j-1):k*j+1,(T0-p)*(t-1)+1:(T0-p)*t)=TS(T0*(t-1)+p-j+1:T0*t-j,:)';
    end
end

% Filling up Y from the data TS
% Filling up Z
% First row (intercept)
% Other rows
% Identifying AR parameters
A   = (Y*Z')/(Z*Z');
% Computing residuals
E   = Y-A*Z;
G=inv(Z*Z');
%estimation of the covariance matrix of the white noise
covE=(Y*Y'-Y*Z'/(Z*Z')*Z*Y')/(size(Y,2)-k*p-1);
if ~exist('Netlabel')  ||isempty(Netlabel) || length(unique(Netlabel))==k 
for i=1:k
    for j=1:k
        index=j+1:k:j+1+(p-1)*k;
        %compute GC,which is equal to the F statistics in conditional
        %Granger Causality analysis
        GCstrength(i,j)=A(i,index)/(G(index,index))*A(i,index)'/covE(i,i);
    end
end
pF=1-fcdf(GCstrength,p,size(Y,2)-k*p-1);%pvalue of F test
pWald=1-chi2cdf(GCstrength,p);%pvalue of Wald test
alpha=0.05;
GC_nonzero=pF<alpha;% We use Bonferroni correction
GCstrength=GCstrength/p;%normalized to granger causality strenth
elseif length(unique(Netlabel))~=k
    for i=1:k diagE(i,i)=covE(i,i);end
    H=kron(diagE,G);
    N=max(Netlabel);
    Wstatistic=zeros(N);
    GCstrength=zeros(N);
    pF=zeros(N);
    pWald=zeros(N);
    for j1=1:N
        for j2=1:N
            roifrom=find(Netlabel==j2);
            roito=find(Netlabel==j1);
            B=zeros(k);B(roito,roifrom)=1;
            for i=1:k B(i,i)=0;end
            b=zeros(k,1);
            for i=1:p b=[b,B];end
            index=find(reshape(b',[],1)==1);
            Aline=reshape(A',[],1);
            if length(index)>0 
                Wstatistic(j1,j2)=Aline(index)'/H(index,index)*Aline(index);%compute the Wald statictics
                GCstrength(j1,j2)=Wstatistic(j1,j2)/length(index);%normalized to granger causality strenth
                pWald(j1,j2)=1-chi2cdf(Wstatistic(j1,j2),length(index));%compute the statistical signifance               
                pF(j1,j2)=1-fcdf(Wstatistic(j1,j2),length(index),size(Y,2)-k*p-1);
            end                
        end
    end
    alpha=0.05;
    GC_nonzero=pF<alpha/N/N;
end
end