%% 1 to X CR
ind=find(CR_ind_up(:,1)==1);ind2=CS_related_hab_ind;
inter=intersect(ind,ind2);
%inter=setdiff(inter,ind_cluster_regressor_checked_cut);
figure,
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'r','filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'b','filled');
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),5,'k','filled');
axis equal;legend('up CR in hab','CS related hab ind','overlap')
figure,
a=[0 0.05];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
plot(squeeze(mean(A_r(:,[trial.hab(2):trial.hab(2)+5],inter),3)))%trial.test(2):trial.test(3)
figure,imagesc(A(:,inter)',[0 0.05])
%clustering
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
M=A_r(frame_ind,:,inter);M=reshape(M,length(frame_ind)*trial.total,[])';
M_norm=normalize(M,2,'zscore');M=M_norm;
%%%clustering
numK=[];[gIXx,~,numK,~]=find_best_k_in_range(M,2:50);
cIXx=1:length(inter);cIX_reg = (1:size(M,1))';
cIX_masterthres={};gIX_masterthres={};kk=1;idd={};
for masterthres=0.7
    para=struct;para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);
    para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
    for ii=1
        [cIX_1toX,gIX_1toX] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
        length(unique(gIX_1toX))
        cIX_reg = (1:size(M,1))';
    end
    length(unique(gIX_1toX))
%     cIX_masterthres{kk}=cIX_1toX;
%     gIX_masterthres{kk}=gIX_1toX;
%     [h]=pushbutton_popupplot_Callback(M,cIX_1toX,gIX_1toX,GetColormap('hsv_new',length(length(unique(gIX_1toX)))),[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
%     title(num2str(masterthres))
%     for zz=unique(gIX_1toX')
%         idd{kk}{zz}=inter(cIX_1toX(find(gIX_1toX==zz)));
%         figure,sepplot(1:size(A,1),A(:,idd{kk}{zz}),GetColormap('hsv_new',length(idd{kk}{zz})));hold on;
%         y=[-12 2];
%         patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
%             repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
%             colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
%         line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))],...
%             [min(y) max(y)],'color','r','linestyle','--','linewidth',1);hold on;
%         xlim([1 size(A,1)])
%          title(['masterthres ' num2str(masterthres) ' cluster ' num2str(zz)])
%     end
%     kk=kk+1;
end

% for kk=1:5
%      cIX_1toX=cIX_masterthres{kk};
%      gIX_1toX=gIX_masterthres{kk};
%     [h]=pushbutton_popupplot_Callback(M,cIX_1toX,gIX_1toX,GetColormap('hsv_new',length(length(unique(gIX_1toX)))),[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
%     title(num2str(masterthres))
%     for zz=unique(gIX_1toX')

%         idd{kk}{zz}=inter(cIX_1toX(find(gIX_1toX==zz)));
%         figure,sepplot(1:size(A,1),A(:,idd{kk}{zz}),GetColormap('hsv_new',length(idd{kk}{zz})));hold on;
%         y=[-12 2];
%         patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
%             [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
%             repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
%             colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
%         line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))],...
%             [min(y) max(y)],'color','r','linestyle','--','linewidth',1);hold on;
%         xlim([1 size(A,1)])
%          title(['masterthres ' num2str(masterthres) ' cluster ' num2str(zz)])
%     end
% end

clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_1toX));
vi='on';
for iii=1
   cIX_iii=cIX_1toX;gIX_iii=gIX_1toX;clrmap_iii=clrmap;
%    cIX_1toX(find(gIX_1toX==1 | gIX_iii==2 | gIX_1toX==3 | gIX_1toX==16 | gIX_1toX==13))=[]
%    gIX_1toX(find(gIX_1toX==1 | gIX_iii==2 | gIX_1toX==3 | gIX_1toX==16 | gIX_1toX==13))=[]
%     [h,gIX2, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M)
%     length(unique(gIX2))
    %2
    [h]=pushbutton_popupplot_Callback(M,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_1toX))
        a(:,ii)=mean(A(:,inter(cIX_1toX(find(gIX_1toX==ii)))),2);
        length(find(gIX_1toX==ii))
    end
    figure,sepplot([1:size(A,1)].*fs.ca,a,GetColormap('hsv_new',max(gIX_1toX)),98);hold on;
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
end

indd_1toX={};indd_1toX_regressor={};
figure,
%scatter3(env.width-supervoxel(ind3,1),supervoxel(ind3,2),supervoxel(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
for ii=1:length(unique(gIX_1toX))
    M_all=A_r(frame_ind,:,:);M_all=reshape(M_all,length(frame_ind)*trial.total,[]);
    [stimcorr,~] = MotorSourceCorrelation(M_all',mean(M(cIX_1toX(find(gIX_1toX==ii)),:),1),[]);
    thr=max(mean(stimcorr)+3.5*std(stimcorr),0.5)
    %figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
    indd_1toX{ii}=intersect(find(stimcorr>thr), find(CR_ind_up(:,1)==1))';
    indd_1toX_regressor{ii}=inter(cIX_1toX(find(gIX_1toX==ii)));
    %figure,plot(a(:,ii));hold on;plot(mean(A(:,indd),2),'r');legend('regressor','stim output')
    subplot(4,4,ii);
    scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),2,[0.5 0.5 0.5],'filled');hold on
    scatter3(supervoxel(indd_1toX{ii},1),supervoxel(indd_1toX{ii},2),supervoxel(indd_1toX{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    title(num2str(ii));%legend;
    view(-90,90)
end

%% cluster间的correlation
a=[];kk=1;
for ii=1:length(unique(gIX_1toX))
    a(:,kk:kk+length(find(gIX_1toX==ii))-1)=A(:,inter(cIX_1toX(find(gIX_1toX==ii))));
    kk=kk+length(find(gIX_1toX==ii));
end
c=corr(a);
figure,imagesc(c,[0 0.5]);colorbar;

%% go-trial 和 nogo trial
go_ind=[37 38 40 41 42];nogo_ind=[39 43 44 45 46];
clr=GetColormap('hsv_new',max(gIX_1toX));
%按CS对齐
figure,subplot(3,1,1)
for ii=1:length(unique(gIX_1toX))
    a=mean(squeeze(mean(A_r(:,go_ind,indd_1toX{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('Go trial');%legend;
subplot(3,1,2)
for ii=1:length(unique(gIX_1toX))
    a=mean(squeeze(mean(A_r(:,nogo_ind,indd_1toX{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('NoGo trial');%legend;
%按运动对齐
subplot(3,1,3)
for ii=1:length(unique(gIX_1toX))
    a=[];kk=1;
    for jj=go_ind
        aa=find(re_startpoint(:,1)==jj & re_startpoint(:,2)>=frameb.cs_start & re_startpoint(:,2)<=frameb.cs_end);
        a(:,kk,:)=A(round(startpoint(aa)/fs.behavior/fs.ca)-frame.per_cycle/2:round(startpoint(aa)/fs.behavior/fs.ca)+frame.per_cycle/2-1,...
           indd_1toX{ii});
       kk=kk+1;
    end
    plot(mean(squeeze(mean(a,3)),2),'color',clr(ii,:),'linewidth',2);hold on;
end
title('Go trial behav');legend;
%% 统计脑区占比
x=[];y=[];
for ii=1:length(unique(gIX_1toX))
    x=[x,indd_1toX{ii}];
    y=[y,ii*ones(length(indd_1toX{ii}),1)'];
end
y=[y,(max(y)+1)*ones(length(setdiff(1:size(A,2),x)),1)']';
x=[x setdiff(1:size(A,2),x)];
clr=GetColormap('hsv_new',length(unique(gIX_1toX))+1);
load('E:\A_Data_lightsheet\Data_huc\20190609\fish3\brain arearegion_mask.mat')
[~,fraction_in_region_in_clust_c,loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],clr,true);
title(['Fish ' num2str(1,'%02d')]);