function ind=indsequence(ind)
global Ncluster;
for i=1:Ncluster
    l(i)=length((find(ind==i)));
    loc(i,1:l(i))=find(ind==i);
end
l2=sort(l);
for i=1:length(l)
    locc=loc(i,:);locc(:,find(locc==0))=[];
    if length(find(l2==l(i)))~=1
        a=find(l2==l(i));
        b=find(l==l2(a(1)));
        for j=1:length(a)            
            locc=loc(b(j),:);locc(:,find(locc==0))=[];
            ind(locc,:)=a(j);
        end
    else
    ind(locc,:)=find(l2==l(i));
    end
end