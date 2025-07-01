function [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,id_in_region_cell,num_in_region_in_clust]=get_region_fraction(region_mask,reg_name,reg_loc,gIX,cIX,loc_all,loc,cmap,isplot)
loc_in_region_in_clust={};loc_in_region_cell={};fraction_in_region_in_clust=nan(size(reg_name,2),2);id_in_region_cell={};
num_in_region_in_clust={};num_in_region_in_clust{1}=nan(size(reg_name,2),2);num_in_region_in_clust{2}=nan(size(reg_name,2),2);
loc_ind=sub2ind(size(region_mask), loc(:,1),loc(:,2),loc(:,3));
% a=region_mask;a(loc_ind(find(loc(:,3)==5)))=30;figure,show_spv_GUI(a)
% figure,scatter(loc(:,1),loc(:,2));axis equal;
% hold on;scatter(reg_loc.R_Hindbrain(2,:),reg_loc.R_Hindbrain(1,:));axis equal;
for rr=1:size(reg_name,2)
    if ~isfield(reg_loc,reg_name{rr})
        loc_in_region_pixel=getfield(reg_loc,strrep(reg_name{rr},' ','_'))';
    else
        loc_in_region_pixel=getfield(reg_loc,reg_name{rr})';
    end
    %figure,scatter3(loc_in_region_pixel(:,1),loc_in_region_pixel(:,2),loc_in_region_pixel(:,3));xlim([1 size(region_mask,1)]);ylim([1 size(region_mask,2)]);axis equal;hold on;
    loc_in_region_pixel=sub2ind(size(region_mask), loc_in_region_pixel(:,2),loc_in_region_pixel(:,1),loc_in_region_pixel(:,3));
    %a=region_mask;a(loc_in_region_pixel)=30;figure,show_spv_GUI(a)
    [loc_in_region_cell{rr},~,id_in_region_cell{rr}]=intersect(loc_in_region_pixel,loc_ind);
    %scatter3(loc(id_in_region_cell{rr},1),loc(id_in_region_cell{rr},2),loc(id_in_region_cell{rr},3));hold on;axis equal;
    loc_in_clust_all=[loc_all(cIX,1),loc_all(cIX,2),loc_all(cIX,3)];
    loc_in_clust_all=sub2ind(size(region_mask),loc_in_clust_all(:,1),loc_in_clust_all(:,2),loc_in_clust_all(:,3));
    for ii=unique(gIX)'
        id_in_clust_in_region=cIX(find(gIX==ii));
        loc_in_clust=[loc_all(id_in_clust_in_region,1),loc_all(id_in_clust_in_region,2),loc_all(id_in_clust_in_region,3)];
        
        %         a=[];
        %         for zz=1:length(loc_in_region_pixel(:,2))
        %             a=find(loc_in_clust(zz,1)==loc_in_region_pixel(:,2) &  loc_in_clust(zz,2)==loc_in_region_pixel(:,1) & loc_in_clust(zz,3)==loc_in_region_pixel(:,3));
        %             a=a+length(a);
        %         end
        loc_in_clust=sub2ind(size(region_mask), loc_in_clust(:,1),loc_in_clust(:,2),loc_in_clust(:,3));
        [loc_in_region_in_clust{rr,ii},~,~]=intersect(loc_in_region_cell{rr},loc_in_clust);
        [loc_in_region_all_clust{rr,ii},~,~]=intersect(loc_in_region_cell{rr},loc_in_clust_all);
        fraction_in_region_in_clust(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(loc_in_region_cell{rr});
        %fraction_in_region_in_clust(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(loc_in_region_all_clust{rr,ii});
        num_in_region_in_clust{1}(rr,ii)=length(loc_in_region_in_clust{rr,ii})/length(id_in_clust_in_region);
        num_in_region_in_clust{2}(rr,ii)=length(loc_in_region_cell{rr});
    end
end
% sum(fraction_in_region_in_clust(:,unique(gIX)),2);a=region_mask;
% for rr=1:size(reg_name,2)
%      for ii=unique(gIX)'
%          a(loc_in_region_in_clust{rr,ii})=30;
%          %a=a+length(loc_in_region_in_clust{rr,ii});
%      end
% end
% figure,show_spv_GUI(a)
%%correct reg_nam
for rr=1:size(reg_name,2)
    reg_name{rr}=lower(reg_name{rr});
%     if strcmp(reg_name{rr},'l pallium rl')
%         reg_name{rr}='l pallium rostral';
%     elseif strcmp(reg_name{rr},'r pallium rl')
%         reg_name{rr}='r pallium rostral';
%     elseif strcmp(reg_name{rr},'l palllium middle')
%         reg_name{rr}='l pallium middle';
%     elseif strcmp(reg_name{rr},'r tectum lateral')
%         reg_name{rr}='r tectum caudal';
%     elseif strcmp(reg_name{rr},'l tectum lateral')
%         reg_name{rr}='l tectum caudal';
%     elseif strcmp(reg_name{rr},'l palllium lateral')
%         reg_name{rr}='l pallium lateral';
%     elseif strcmp(reg_name{rr},'logitude')
%         reg_name{rr}='longitude';
%     end
    if strcmp(reg_name{rr},'lcb')
        reg_name{rr}='lCb';
    elseif strcmp(reg_name{rr},'lhb')
        reg_name{rr}='lHb';
    elseif strcmp(reg_name{rr},'lp')
        reg_name{rr}='lPa';
    elseif strcmp(reg_name{rr},'lt')
        reg_name{rr}='lOT';
    elseif strcmp(reg_name{rr},'rcb')
        reg_name{rr}='rCb';
    elseif strcmp(reg_name{rr},'rhb')
        reg_name{rr}='rHb';
    elseif strcmp(reg_name{rr},'rp')
        reg_name{rr}='rPa';
    elseif strcmp(reg_name{rr},'rt')
        reg_name{rr}='rOT';
            elseif strcmp(reg_name{rr},'h')
        reg_name{rr}='H';
    elseif strcmp(reg_name{rr},'tl')
        reg_name{rr}='TL';
    end

end
% Lable={'l pallium rostral'; 'r pallium rostral';'l pallium middle';'r pallium middle';'l pallium lateral';'r pallium lateral';...
%     'l habenula';'r habenula'; 'longitude';...
%     'l tectum rostral';'r tectum rostral';'l tectum middle'; 'r tectum middle';'l tectum caudal';'r tectum caudal';...
%     'l tegmentum';'r tegmentum';...
%     'l cerebellum';'r cerebellum';...
%     'l hindbrain';'r hindbrain';...
%    
%     };
Lable={'rPa';'lPa'; 'rHb';'lHb';'TL';'rOT';'lOT';'rCb';'lCb';'H';};
id=[];
for rr=1:size(Lable,1)
    id(rr)=find(strcmp(Lable{rr},reg_name));
end
type=2;
if isplot
    figure('position',[200,200,1300,800]),
    if ~isempty(fraction_in_region_in_clust)
        c = categorical({reg_name{id}})';
        switch type
            case 1
                %id=[6,17,5,16,4,15,2,13,9,20,8,19,7,18,10,21,1,12,3,14,11];%[1:size(reg_name,2)];%[1,8,2,9,3,10,4,11,5,12,6,13,7];
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

fraction_in_region_in_clust=fraction_in_region_in_clust(id,:);
num_in_region_in_clust{1}=num_in_region_in_clust{1}(id,:);
num_in_region_in_clust{2}=num_in_region_in_clust{2}(id,:);

