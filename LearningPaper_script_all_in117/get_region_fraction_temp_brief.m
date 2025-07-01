function [loc_in_region_in_clust_r,index_in_region_in_clust_r,fraction_in_region_in_clust_r,brain_region_id,loc_in_region_in_temp,id_in_region_in_temp,num_in_region_in_clust_r,Lable]=get_region_fraction_temp_brief(gIX,ind_supervoxel,supervoxel,loc_temp,region_mask,reg_name,cmap,isplot)
index_in_region_in_clust={};brain_region_id=strings(length(supervoxel),2);
loc_in_region_in_clust={};loc_in_region_in_temp={};fraction_in_region_in_clust=nan(size(reg_name,2),max(unique(gIX)));id_in_region_in_temp={};
num_in_region_in_clust={};num_in_region_in_clust{1}=nan(size(reg_name,2),max(unique(gIX)));num_in_region_in_clust{2}=nan(size(reg_name,2),max(unique(gIX)));num_in_region_in_clust{3}=nan(size(reg_name,2),max(unique(gIX)));

%figure,imagesc(region_mask(:,:,20))

if ~isempty(supervoxel)
    ind=find(supervoxel(:,1)<=0 | supervoxel(:,2)<=0 | supervoxel(:,3)<=0 |supervoxel(:,1)>size(region_mask,2) | supervoxel(:,2)>size(region_mask,1) | supervoxel(:,3)>size(region_mask,3));
    gIX(ind)=[];
    supervoxel(ind,:)=[];
    ind_supervoxel(ind)=[];
end
ind_temp=1:size(loc_temp,1);
if ~isempty(loc_temp)
    ind=find(loc_temp(:,1)<=0 | loc_temp(:,2)<=0 | loc_temp(:,3)<=0 |loc_temp(:,1)>size(region_mask,2) | loc_temp(:,2)>size(region_mask,1) | loc_temp(:,3)>size(region_mask,3)-1);
    loc_temp(ind,:)=[];
    ind_temp(ind)=[];
end
redion_id=nan(size(loc_temp,1),1);
for ii=1:size(loc_temp,1)
redion_id(ii)=region_mask(round(loc_temp(ii,2)),round(loc_temp(ii,1)),ceil(loc_temp(ii,3)));
end
for rr=1:size(reg_name,2)
    id_in_region_in_temp{rr}=ind_temp(find(redion_id==rr));
    loc_in_region_in_temp{rr}=loc_temp(find(redion_id==rr),:);
end
%figure,scatter(supervoxel(:,1),supervoxel(:,2));hold on;scatter(supervoxel(id_in_region_cell{rr},1),supervoxel(id_in_region_cell{rr},2));hold on
for ii=unique(gIX)'
    id_in_clust_in_region=find(gIX==ii);
    loc_in_clust=[supervoxel(id_in_clust_in_region,1),supervoxel(id_in_clust_in_region,2),supervoxel(id_in_clust_in_region,3)];
    redion_id=nan(size(id_in_clust_in_region,2),1);
    for zz=1:length(id_in_clust_in_region)
        redion_id(zz)=region_mask(round(loc_in_clust(zz,2)),round(loc_in_clust(zz,1)),ceil(loc_in_clust(zz,3)));
    end  
    for rr=1:size(reg_name,2)
        index_in_region_in_clust{rr,ii}=ind_supervoxel(id_in_clust_in_region(find(redion_id==rr)));
        loc_in_region_in_clust{rr,ii}=loc_in_clust(find(redion_id==rr),:);
        fraction_in_region_in_clust(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(loc_in_region_in_temp{rr});
        num_in_region_in_clust{1}(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(id_in_clust_in_region);
        num_in_region_in_clust{2}(rr,ii)=length(loc_in_region_in_clust{rr,ii});
        num_in_region_in_clust{3}(rr,ii)=length(id_in_clust_in_region);
        brain_region_id(index_in_region_in_clust{rr,ii},1)=rr;
        brain_region_id(index_in_region_in_clust{rr,ii},2)=string(reg_name{rr});
    end
                
end


Lable={'L OB','R OB','L P','R P',...
    'L H','R H','L TeO','R TeO','L Np','R Np','L Otr','R Otr','TL','PO',...
    'Pr','PT','Th','rH','iH','cH','mPT','T',...
    'L TS','R TS','Va','CCe','L gT','R gT','R1','R2'};

id=[];
for rr=1:size(Lable,2)
    id(rr)=find(strcmp(strrep(Lable{rr},'_',''),reg_name));
end
type=1;
if isplot
    figure('position',[200,200,1300,800]),
    if ~isempty(fraction_in_region_in_clust)
        c = categorical({reg_name{id}})';
        switch type
            case 1
                y=[fraction_in_region_in_clust(id,:)];
                b=bar(y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]');
                %b=bar(c',y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]');
                set(gca,'xticklabels',{reg_name{id}},'XTickLabelRotation',45,'fontsize',12);
                %cmap(end+1:max(unique(gIX)),:)= cmap(end,:);
                kk=1;
                for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
                    b(k).CData = cmap(kk,:);
                    kk=kk+1;
                end
                legend([b(unique(gIX))],num2str(unique(gIX)),'location','bestoutside');ylim([0 1]);%
            case 2
                y=num_in_region_in_clust{2}(id,1);
                b=bar(y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);hold on;
                %b=bar(c,y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);hold on;
                b.CData = [0.5 0.5 0.5];
                y=num_in_region_in_clust{1}(id,:);
                b=bar(y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);hold on;
                %b=bar(c,y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);hold on;
                set(gca,'xticklabels',{reg_name{id}},'XTickLabelRotation',45,'fontsize',12);
                kk=1;
                for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
                    b(k).CData = cmap(kk,:);
                    kk=kk+1;
                end               
                legend([b(unique(gIX))],num2str(unique(gIX)),'location','bestoutside');
                text([1:length(c)],num_in_region_in_clust{2}(id,1)+0.1*num_in_region_in_clust{2}(id,1),num2str(fraction_in_region_in_clust(id,unique(gIX)),'%.2f\n'));
        end
    else
        disp('is empty!!')
    end
end

fraction_in_region_in_clust_r=fraction_in_region_in_clust(id,:);
num_in_region_in_clust_r{1}=num_in_region_in_clust{1}(id,:);
num_in_region_in_clust_r{2}=num_in_region_in_clust{2}(id,:);
num_in_region_in_clust_r{3}=num_in_region_in_clust{3}(id,:);
index_in_region_in_clust_r={};loc_in_region_in_clust_r={};
for ii=1:size(index_in_region_in_clust,2)
    for jj=1:length(id)
    index_in_region_in_clust_r{jj,ii}=index_in_region_in_clust{id(jj),ii};
    loc_in_region_in_clust_r{jj,ii}=loc_in_region_in_clust{id(jj),ii};
    end
end

