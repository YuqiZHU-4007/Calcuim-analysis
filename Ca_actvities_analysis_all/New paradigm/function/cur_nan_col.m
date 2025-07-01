function X=cur_nan_col(X)

kk=1;cut=[];
for ii=1:size(X,1)
    if sum(isnan(X(ii,:)))==size(X,2)
       cut(kk)=ii;kk=kk+1;
    end
end
X(cut,:)=[];