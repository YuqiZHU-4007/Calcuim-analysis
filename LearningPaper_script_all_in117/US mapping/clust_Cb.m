%% clust Cb
load('E:\A_Data_lightsheet\Data_huc\20190609\fish3\brain arearegion_mask.mat')
y=[ones(size(A,2),1)]';
x=[1:size(A,2)];
[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,id_in_region_cell,num_in_region_in_clust]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],[0.5 0.5 0.5],true);
title(['Fish ' num2str(1,'%02d')]);

reg_ind=[1,12];reg_id=[];
for ii=1:length(reg_ind)
    reg_id=[reg_id;id_in_region_cell{reg_ind(ii)}];
end
ind=find(CR_ind_up(:,1)==1);ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);ind2=intersect(ind2, CS_related_hab_down_ind);
ind3=union(ind,ind2);ind3=setdiff(1:size(A,2),ind3);
ind=find(CR_ind_up(:,2)==1);
ind=intersect(ind,CS_related_tst_ind);
inter=intersect(ind,ind3);
reg_id=intersect(inter,reg_id);
figure,imagesc(A(:,reg_id)',[0 0.05])

figure,
scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),2,[0.5 0.5 0.5],'filled');hold on
scatter3(supervoxel(reg_id,1),supervoxel(reg_id,2),supervoxel(reg_id,3))

%clustering
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
M=A_r(frame_ind,:,reg_id);M=reshape(M,length(frame_ind)*trial.total,[])';
M_norm=normalize(M,2,'zscore');M=M_norm;
%cluster
numK=[];[gIXx,~,numK,~]=find_best_k_in_range(M,5:50);
cIXx=1:length(reg_id);cIX_reg = (1:size(M,1))';
masterthres=0.6;
para=struct;para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',5);
for ii=1
[cIX_Cb,gIX_Cb] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX_Cb))
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_Cb));
vi='on';
for iii=1
    cIX_iii=cIX_Cb;gIX_iii=gIX_Cb;clrmap_iii=clrmap;
    B=M;
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_Cb))
        a(:,ii)=mean(A(:,reg_id(cIX_Cb(find(gIX_Cb==ii)))),2);
        length(find(gIX_Cb==ii))
    end
    figure,sepplot([1:size(A,1)].*fs.ca,a,GetColormap('hsv_new',max(gIX_Cb)));hold on;
    y=[-12 2];
    patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']'.*fs.ca,...
        repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
     line([startpoint;startpoint]./fs.behavior,...
        [min(y) max(y)],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',0.5);hold on;    
    line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))].*fs.ca,...
        [min(y) max(y)],'color','r','linestyle','--','linewidth',2);hold on;
    xlim([1 size(A,1)].*fs.ca)  
    
    figure,
    scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),2,[0.5 0.5 0.5],'filled');hold on
    for ii=1:length(unique(gIX_Cb))
        a=find(gIX_Cb==ii);
        scatter3(supervoxel(reg_id(cIX_Cb(a)),1),supervoxel(reg_id(cIX_Cb(a)),2),supervoxel(reg_id(cIX_Cb(a)),3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    axis equal;
end