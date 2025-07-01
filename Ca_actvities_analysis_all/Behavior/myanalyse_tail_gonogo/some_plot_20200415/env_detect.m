function [env_num,env_loc_m,env_end_m,m_x]=env_detect(x,bin,type)
env_num=0;
env_loc=[];
env_end=[];
switch type
    case 1
        a=find(x~=0);
        if ~isempty(a)
        dif=diff(a);
        ind=find(dif~=1)+1;%ind=find(dif~=1 & dif>=bin)+1;        
        env_num=length(ind)+1;
        env_loc=max([max(a(1),1);a(ind)]-1,1);        
        ind=find(dif~=1);%ind=find(dif~=1 & dif>=bin);
        env_end=[a(ind);a(end)];
        end
        %cut baseline
        for ii=1:env_num
            ind=env_loc(ii):env_end(ii);
            a=x(ind);aa=abs(diff([a(1);a]));
            ind_cut=find(aa<0.5 & [aa(2:end);aa(end)]<0.5);
            ind(ind_cut)=[];
            if ~isempty(ind)
                env_loc(ii)=min(ind);
                a=x(ind);aa=abs(diff([a;a(end)]));
                if ~isempty(aa)
                    ind_cut=find([aa(1);aa(1:end-1)]<0.5 & aa <0.5);
                    ind(ind_cut)=[];
                    env_end(ii)=max(ind);
                end
            else
                env_loc(ii)=nan;
                env_end(ii)=nan;
            end
        end
    case 2
        dif=diff(x);
end
ind=isnan(env_loc);
env_loc(isnan(env_loc))=[];
env_end(isnan(env_end))=[];
%merge 
env_loc_m=env_loc;env_end_m=env_end;
for ii=1:length(env_loc)-1
    for jj=ii+1:length(env_loc)
        if (env_loc(jj)-env_end_m(ii))<bin
            env_end_m(ii)=env_end(jj);env_end(jj);
            env_loc_m(jj)=nan;
        end
    end
end
% ind=find((env_loc(2:end)-env_end(1:end-1))<bin);
% env_loc_m(ind(ii)+1)=[];
% env_end_m(ind)=env_end(ind+1);env_end_m(ind+1)=[];
env_end_m(isnan(env_loc_m))=[];env_loc_m(isnan(env_loc_m))=[];
env_dur=env_end_m-env_loc_m; env_end_m(find(env_dur<=1))=[];env_loc_m(find(env_dur<=1))=[];
env_num=length(env_loc_m);
m_x=[];
for ii=1:length(env_loc_m)
    m_x(ii)=mean(x(env_loc_m(ii):env_end_m(ii)));
end
% figure,plot(x);hold on;scatter(env_loc,x(env_loc),'r','filled');hold on;
% scatter(env_end,x(env_end),'k','filled');hold on;
% scatter(env_loc_m,x(env_loc_m),'y','filled');hold on;
% scatter(env_end_m,x(env_end_m),'g','filled');hold on;