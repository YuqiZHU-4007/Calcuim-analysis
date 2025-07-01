function tind=detect_tail(text,ind,indend,base2,base3)%%%%text 为taildeg_everytrial（i，：）
%%%tind 为duration
base2=0.4*10;
l=length(ind);
tind=zeros;
if l==1
    dif=diff(text(ind:indend));
    tind=max(find([dif(2:end) dif(end)]<base2 & [dif(1) dif(1:end-1)]>=base2 & dif>=base2 ));%& [dif(3:end) dif(end) dif(end)]<base2
    if length(tind)~=0
        tind=tind+1+ind-ind;
    else
        tind=indend-ind+1;
    end
    if ind+tind+20<=indend
        if length(find(text(ind+tind-6:ind+tind+20)==0))~=0
            if length(find(text(ind+tind-6:ind+tind)))~=0
                h=find(text(ind+tind-6:ind+tind)==0);
                tind=tind-h(1);
            else
                h=find(text(ind+tind:ind+tind+20)==0);
                tind=tind+h(1);
            end
        end
    else
        if length(find(text(ind+tind-6:indend)==0))~=0
            if length(find(text(ind+tind-6:indend)))~=0
                h=find(text(ind+tind-6:ind+tind)==0);
                tind=tind-h(1);
            else
                h=find(text(ind+tind:indend)==0);
                tind=tind+h(1);
            end
        end
    end
else
    for z=1:l
        if z==l
            dif=diff(text(ind(z):indend(end)));
        else
            dif=diff(text(ind(z):ind(z+1)-10));
        end
        if size(dif,1)<=5
            a=size(dif,1);
        else
            a=max(find([dif(2:end) dif(end)]<base2 & [dif(1) dif(1:end-1)]>=base2 & dif>=base2 & [dif(3:end) dif(end) dif(end)]<base2 & [dif(4:end) dif(end) dif(end) dif(end)]<base2));
        end
        %taildeg_everytrial(i,mind(z)-bin2:mind(z)-1)=0;
        if length(a)~=0
            tind(z)=a;
            tind(z)=max(tind(z)+1+ind(z)-ind(z),min(indend(z)+20,indend(end))+1-ind(z));%%%25
        else
            if z==l
                tind(z)=indend(end)-ind(z)+1;
            else
                %tind(z)=ind(z+1)-5-ind(z)+1;%%%判断不出尾
                tind(z)=min(indend-base3*(l-z)-ind(z)+1,ind(z+1)-base3-ind(z)+1,indend(z)+25+1-ind(z));
            end
        end
        if ind(z)+tind(z)+20<=indend(end)
            if length(find(text(ind(z)+tind(z)-6:ind(z)+tind(z)+20)==0))~=0
                if length(find(text(ind(z)+tind(z)-6:ind(z)+tind(z))))~=0
                    h=find(text(ind(z)+tind(z)-6:ind(z)+tind(z))==0);
                    tind(z)=tind(z)-h(1);
                else
                    h=find(text(ind(z)+tind(z):ind(z)+tind(z)+20)==0);
                    tind(z)=tind(z)+h(1);
                end
            end
        else
            if length(find(text(ind(z)+tind(z)-6:indend(end))==0))~=0
                if length(find(text(ind(z)+tind(z)-6:indend(end))))~=0
                    h=find(text(ind(z)+tind(z)-6:ind(z)+tind(z))==0);%%%%%
                    tind(z)=tind(z)-h(1);
                else
                    h=find(text(ind(z)+tind(z):indend(end))==0);
                    tind(z)=tind(z)+h(1);
                end
            end
        end
    end
end
