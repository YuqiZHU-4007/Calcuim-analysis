function onset=findonset(A,thr)
n=3; onset=nan(1,size(A,2));
for ii=1:size(A,2)
    for t=1:size(A,1)-n
        if A(t,ii)>thr && sum(A(t:t+n,ii)>thr)>=n-1
            onset(:,ii)=t;
        else
            break;
        end
    end
end
onset(find(onset==0))=nan;
end
