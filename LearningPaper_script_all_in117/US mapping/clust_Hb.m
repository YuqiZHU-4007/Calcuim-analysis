%% clust Hb
load('E:\A_Data_lightsheet\Data_huc\20190609\fish3\brain arearegion_mask.mat')
y=[ones(size(A,2),1)]';
x=[1:size(A,2)];
[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,id_in_region_cell,num_in_region_in_clust]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],[0.5 0.5 0.5],true);
title(['Fish ' num2str(1,'%02d')]);

reg_ind=[2];reg_id=[];
for ii=1:length(reg_ind)
    reg_id=[reg_id;id_in_region_cell{reg_ind(ii)}];
end
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
masterthres=0.7;
para=struct;para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);
para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);
para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
for ii=1
[cIX_Hb,gIX_Hb] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX_Hb))
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_Hb));
vi='on';
for iii=1
    cIX_iii=cIX_Hb;gIX_iii=gIX_Hb;clrmap_iii=clrmap;
    B=M;
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_Hb))
        a(:,ii)=mean(A(:,reg_id(cIX_Hb(find(gIX_Hb==ii)))),2);
        length(find(gIX_Hb==ii))
    end
    figure,sepplot([1:size(A,1)].*fs.ca,a,GetColormap('hsv_new',max(gIX_Hb)));hold on;
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
    for ii=1:length(unique(gIX_Hb))
        a=find(gIX_Hb==ii);
        scatter3(supervoxel(reg_id(cIX_Hb(a)),1),supervoxel(reg_id(cIX_Hb(a)),2),supervoxel(reg_id(cIX_Hb(a)),3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    axis equal;
end

%% cluster间的correlation
a=[];kk=1;
for ii=1:length(unique(gIX_Hb))
    a(:,kk:kk+length(find(gIX_Hb==ii))-1)=A(:,reg_id(cIX_Hb(find(gIX_Hb==ii))));
    kk=kk+length(find(gIX_Hb==ii));
end
c=corr(a);
figure,imagesc(c,[0 0.5]);colorbar;

%% go-trial 和 nogo trial
go_ind=[37 38 40 41 42];nogo_ind=[39 43 44 45 46];
clr=GetColormap('hsv_new',max(gIX_0to1));
%按CS对齐
figure,subplot(3,1,1)
for ii=1:length(unique(gIX_Hb))
    a=mean(squeeze(mean(A_r(:,go_ind,indd_0to1{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('Go trial');%legend;
subplot(3,1,2)
for ii=1:length(unique(gIX_Hb))
    a=mean(squeeze(mean(A_r(:,nogo_ind,indd_0to1{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('NoGo trial');%legend;
%按运动对齐
subplot(3,1,3)
for ii=1:length(unique(gIX_Hb))
    a=[];kk=1;
    for jj=go_ind
        aa=find(re_startpoint(:,1)==jj & re_startpoint(:,2)>=frameb.cs_start-10 & re_startpoint(:,2)<=frameb.cs_end);
        a(:,kk,:)=A(round(startpoint(aa)/fs.behavior/fs.ca)-frame.per_cycle/2:round(startpoint(aa)/fs.behavior/fs.ca)+frame.per_cycle/2-1,...
           indd_0to1{ii});
       kk=kk+1;
    end
    plot(mean(squeeze(mean(a,3)),2),'color',clr(ii,:),'linewidth',2);hold on;
end
title('Go trial behav');legend;