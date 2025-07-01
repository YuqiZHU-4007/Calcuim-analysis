function [loc_in_region_in_clust,index_in_region_in_clust_r,fraction_in_region_in_clust_r,loc_in_region_cell,id_in_region_cell,num_in_region_in_clust_r,Lable]=get_region_fraction_temp(region_mask,reg_name,reg_loc,gIX,surpervolxel,loc_temp,cmap,isplot,nn)
if ~isempty(strfind(nn,'2019'))
    res2=[0.66,0.66,8];
elseif ~isempty(strfind(nn,'2021')) | ~isempty(strfind(nn,'2022'))
    res2=[0.66,0.66,10];
end
res_temp=[0.66 0.66 10];
index_in_region_in_clust={};
loc_in_region_in_clust={};loc_in_region_cell={};fraction_in_region_in_clust=nan(size(reg_name,2),10);id_in_region_cell={};
num_in_region_in_clust={};num_in_region_in_clust{1}=nan(size(reg_name,2),10);num_in_region_in_clust{2}=nan(size(reg_name,2),10);num_in_region_in_clust{3}=nan(size(reg_name,2),10);

if ~isempty(loc_temp)
    loc_temp(find(loc_temp(:,1)<=0 | loc_temp(:,1)>size(region_mask,1)*res2(1)),:)=[];
end
if ~isempty(loc_temp)
    loc_temp(find(loc_temp(:,2)<=0 | loc_temp(:,2)>size(region_mask,2)*res2(2)),:)=[];
end
if ~isempty(loc_temp)
    loc_temp(find(loc_temp(:,3)<=0 | loc_temp(:,3)>(size(region_mask,3)-1)*res2(3)),:)=[];
end
loc_ind=sub2ind(size(region_mask),floor(loc_temp(:,2)/res_temp(2)),floor(loc_temp(:,1)/res_temp(1)),floor(loc_temp(:,3)/res_temp(3))+1);
% a=region_mask;
% for ii=1:length(loc_ind)
% a(loc_ind(ii))=100;
% end
% figure,show_spv_GUI(a)
% loc_in_clust_all=surpervolxel;
% loc_in_clust_all=sub2ind(size(region_mask),floor(loc_in_clust_all(:,1)/res(1)),floor(loc_in_clust_all(:,2)/res(2)),floor(loc_in_clust_all(:,3)/res(3))+1);
for rr=1:size(reg_name,2)
    if ~isfield(reg_loc,reg_name{rr})
        loc_in_region_pixel=getfield(reg_loc,strrep(reg_name{rr},' ','_'));
    else
        loc_in_region_pixel=getfield(reg_loc,reg_name{rr});
    end
    %figure,scatter3(loc_in_region_pixel(:,1),loc_in_region_pixel(:,2),loc_in_region_pixel(:,3));xlim([1 size(region_mask,1)]);ylim([1 size(region_mask,2)]);axis equal;hold on;
    loc_in_region_pixel_subind=sub2ind(size(region_mask), loc_in_region_pixel(:,2)/res_temp(2),loc_in_region_pixel(:,1)/res_temp(1),loc_in_region_pixel(:,3)/res_temp(3)+1);
    %a=region_mask;a(loc_in_region_pixel_subind)=0;show_spv_GUI(a)
    [loc_in_region_cell{rr},~,id_in_region_cell{rr}]=intersect(loc_in_region_pixel_subind,loc_ind);
    % [loc_in_region_all_clust{rr},~,~]=intersect(loc_in_region_cell{rr},loc_in_clust_all);
    %figure,scatter3(loc_temp(:,1),loc_temp(:,2),loc_temp(:,3));hold on;axis equal;
    %figure,scatter3(loc_temp(id_in_region_cell{rr},1),loc_temp(id_in_region_cell{rr},2),loc_temp(id_in_region_cell{rr},3),12,'r');hold on;axis equal;
    
    for ii=unique(gIX)'
        id_in_clust_in_region=(find(gIX==ii));
        loc_in_clust=[surpervolxel(id_in_clust_in_region,1),surpervolxel(id_in_clust_in_region,2),surpervolxel(id_in_clust_in_region,3)];
%         figure,scatter3(loc_in_clust(:,1),loc_in_clust(:,2),loc_in_clust(:,3));hold on;
        loc_in_clust_subind=sub2ind(size(region_mask), max(floor(loc_in_clust(:,2)/res2(2)),1),max(floor(loc_in_clust(:,1)/res2(1)),1),floor(loc_in_clust(:,3)/res2(3))+1);
%         a=region_mask;b=[floor(loc_in_clust(:,2)/res2(2)),floor(loc_in_clust(:,1)/res2(1)),floor(loc_in_clust(:,3)/res2(3))+1];
%         for zz=1:length(loc_in_clust)
%             a(b(zz,:))=0;sub2ind(size(region_mask), floor(loc_in_clust(zz,2)/res2(2)),floor(loc_in_clust(zz,1)/res2(1)),floor(loc_in_clust(zz,3)/res2(3))+1);
%         end
%         a=region_mask;
%         for zz=1:length(loc_in_clust_subind)
%             a(loc_in_clust_subind(zz))=100;
%         end
%         show_spv_GUI(a)
        [loc_in_region_in_clust{rr,ii},~,ib]=intersect(loc_in_region_pixel_subind,loc_in_clust_subind);
        index_in_region_in_clust{rr,ii}=id_in_clust_in_region(ib);
%         figure,scatter(loc_in_clust(:,1),loc_in_clust(:,2));hold on;
%                 a=region_mask;
%         for zz=1:length(loc_in_region_in_clust{rr,ii})
%             a(loc_in_region_in_clust{rr,ii}(zz))=100;
%         end
%         show_spv_GUI(a)
        fraction_in_region_in_clust(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(loc_in_region_cell{rr});
        %fraction_in_region_in_clust(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(loc_in_region_all_clust{rr,ii});
        num_in_region_in_clust{1}(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(id_in_clust_in_region);
        num_in_region_in_clust{2}(rr,ii)=length(loc_in_region_in_clust{rr,ii});
        num_in_region_in_clust{3}(rr,ii)=length(loc_in_region_cell{rr});
    end
end
% sum(fraction_in_region_in_clust(:,unique(gIX)),2);a=region_mask;

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
index_in_region_in_clust_r={};
for ii=1:size(index_in_region_in_clust,2)
    for jj=1:length(id)
    index_in_region_in_clust_r{jj,ii}=index_in_region_in_clust{id(jj),ii};
    end
end

