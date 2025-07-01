function t_dur=find_t_dur(taildeg_everytrial,trailnum)
t_dur=zeros;
for i=1:trailnum
    indsta=zeros;indend=zeros;
    a=taildeg_everytrial(i,:);
    ind=find(a~=0);
    %     tail_sti=a(sti_start+1:sti_start+duration);
    %     tail_base=a(1:sti_start);
    if length(ind)>1      
        dif=diff(ind);
        uni=unique(dif);
        if length(uni)~=1
            indsta(1)=ind(1);
            k=2;
            for j=2:length(uni)
                if length(find(dif==uni(j)))==1
                    indsta(k)=ind(find(dif==uni(j))+1);
                    indend(k-1)=ind(find(dif==uni(j)));
                    k=k+1;
                else
                    indsta(k:k+length(find(dif==uni(j)))-1)=ind(find(dif==uni(j))+1);
                    indend(k-1:k+length(find(dif==uni(j)))-2)=ind(find(dif==uni(j)));
                    k=k+length(find(dif==uni(j)));
                end
            end
            indend(end+1)=ind(end);
        else
            indsta=ind(1);indend=ind(end);
        end
        indsta=sort(indsta);indend=sort(indend);
        dur=indend-indsta+1;
        t_dur(1:length(indsta),1,i)=indsta';
        t_dur(1:length(indsta),2,i)=dur';
    else
        t_dur(:,:,i)=0;
    end 
end
for i=1:trailnum
    tt=zeros;ind=zeros;
    t=t_dur(:,:,i);t(find(t(:,1)==0),:)=[];
    if size(t,1)~=0
    tt=t+[t(:,2)-1 zeros(size(t,1),1)];tt=[t(2:end,:);t(1,:)]-tt;
    ind=find((tt(:,1))==2);
    if ind~=0
    t_dur(ind,:,i)=[t_dur(ind,1,i) t_dur(ind,2,i)+t_dur(ind+1,2,i)+1];
    t_dur(ind+1,:,i)=0;  
    end
    end
    %t_dur(find(t_dur(:,1,i)==0),:,i)=[];
end