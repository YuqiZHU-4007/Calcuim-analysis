function [X,X1]=updata_struct(X,name,areaa,X1)
if ~isempty(areaa)
    m=mean(areaa,2);
    sd=std(areaa,[],2);
    n=size(areaa,2);
    X=setfield(X,name,[m,sd,n*ones(size(m))]);
    X1=setfield(X1,name,areaa);
else
    X=setfield(X,name,[]);
    X1=setfield(X1,name,[]);
end