%估算脑区null 分布
function [P_z,H_z,P_t,H_t,fraction_samlpe]=statistic_P_in_region_fraction(loc,fraction,all_loc,temp_env,temp_supervoxel,reg_mask,reg_name,reg_loc)
%fraction
%loc
clr=jet(2);N=10000;
sample_num=size(loc.ind,1);
fraction_samlpe=[];fraction_samlpe_parfor={};
%% 估计null distribution
parfor nsample=1:N
    kk=1;
    rand_id=randi(size(all_loc.ind,1),sample_num,1);
    random_n=all_loc(rand_id,:);
    for nn=string(unique(random_n.label))'
        ind=find(strcmp(string(random_n.label),nn));
        surpervolxel=[random_n.x(ind) random_n.y(ind) random_n.z(ind)];
        if ~isempty(strfind(nn,'2019'))
            res=[0.66,0.66,8];
        elseif ~isempty(strfind(nn,'2021'))
            res=[0.66,0.66,10];
        end
        gIX=[ones(size( surpervolxel,1),1)];
        gIX(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res(1)))=[];surpervolxel(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res(1)),:)=[];
        gIX(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res(2)))=[];surpervolxel(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res(2)),:)=[];
        gIX(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res(3)))=[];surpervolxel(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res(3)),:)=[];
        [~,fraction_in_region_in_clust,~,~,~]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
            gIX,surpervolxel,temp_supervoxel,clr,false,nn);
        %fraction_samlpe(:,nsample,kk)=fraction_in_region_in_clust(:,1);kk=kk+1;
        fraction_samlpe_parfor{nsample}(:,kk)=fraction_in_region_in_clust(:,1);kk=kk+1;
    end
end
for ii=1:N
    fraction_samlpe(:,ii,:)=fraction_samlpe_parfor{ii};
end
%% 统计
P_z=[];H_z=[];alpha=0.05;%双尾
for region_i=1:size(fraction,1)
    a=mean(squeeze(fraction_samlpe(region_i,:,:)),2,'omitnan');
    m=mean(a);s=std(a);n=length(fraction(region_i,:));
    t=(mean(fraction(region_i,:))-m)./(s*sqrt(n));
    t_a=norminv(alpha/2);
    if abs(t)>abs(t_a) h=true;else h=false;end
    p = 2 * normcdf(-abs(t));
    P_z(region_i)=p;
    H_z(region_i)=h;
end
% find(P_z<alpha)
% find(H_z==1)

P_t=[];H_t=[];alpha=0.05;%双尾
for region_i=1:size(fraction,1)
    a=mean(squeeze(fraction_samlpe(region_i,:,:)),2,'omitnan');
    m=mean(a);s=std(a);n=length(fraction(region_i,:));
    t=(mean(fraction(region_i,:))-m)./(s*sqrt(n));
    t_a=tinv(alpha/2,n-1);
    if abs(t)>abs(t_a) h=true;else h=false;end
    p = 2 * tcdf(-abs(t),n-1);
    P_t(region_i)=p;
    H_t(region_i)=h;
end
% find(P_t<alpha)
% find(H_t==1)