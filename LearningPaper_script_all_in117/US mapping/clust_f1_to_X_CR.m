%% -1 to X CR 
ind2=find(CR_ind_down(:,1)==1);
inter=intersect(ind2,CS_related_hab_down_ind);
figure,
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),10,'r','filled');
axis equal;legend('0-1 CR ')
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
%cluster
numK=[];[gIXx,~,numK,~]=find_best_k_in_range(M,2:30);
cIXx=1:length(inter);cIX_reg = (1:size(M,1))';
masterthres=0.6;
para=struct;para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',5);
for ii=1
[cIX_f1toX,gIX_f1toX] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX_f1toX))
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_f1toX));
vi='on';
for iii=1
    cIX_iii=cIX_f1toX;gIX_iii=gIX_f1toX;clrmap_iii=clrmap;
    %     [h,gIX2, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M)
    %     length(unique(gIX2))
    %2
    B=M;
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_f1toX))
        a(:,ii)=mean(A(:,inter(cIX_f1toX(find(gIX_f1toX==ii)))),2);
        length(find(gIX_f1toX==ii))
    end
    figure,sepplot(1:size(A,1),a,GetColormap('hsv_new',max(gIX_f1toX)));hold on;
    y=[-12 2];
    patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
        repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))],...
        [min(y) max(y)],'color','r','linestyle','--','linewidth',1);hold on;
    xlim([1 size(A,1)])
    
    figure,
    indd_f1toX={};aa=[];indd_f1toX_regressor={};
    for ii=1:length(unique(gIX_f1toX))
        M_all=A_r(frame_ind,:,:);M_all=reshape(M_all,length(frame_ind)*trial.total,[]);
        [stimcorr,~] = MotorSourceCorrelation(M_all',mean(M(cIX_f1toX(find(gIX_f1toX==ii)),:),1),[]);
         n=3;
        thr=max(mean(stimcorr)+3.5*std(stimcorr),0.5)
        %figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
        indd_f1toX{ii}=intersect(find(stimcorr>thr),find(CR_ind_down(:,1)==1) )';
        indd_f1toX_regressor{ii}=inter(cIX_f1toX(find(gIX_f1toX==ii)));
        aa(:,ii)=mean(A(:,indd_f1toX{ii}),2);
        %figure,plot(a(:,ii));hold on;plot(mean(A(:,indd{ii}),2),'r');legend('regressor','stim output')
        scatter3(supervoxel(indd_f1toX{ii},1),supervoxel(indd_f1toX{ii},2),supervoxel(indd_f1toX{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    figure,
    for ii=1:length(unique(gIX_f1toX))
        scatter3(supervoxel(indd_f1toX_regressor{ii},1),supervoxel(indd_f1toX_regressor{ii},2),supervoxel(indd_f1toX_regressor{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    figure,sepplot(1:size(A,1),aa,GetColormap('hsv_new',max(gIX_f1toX)));hold on;
    y=[-12 2];
    patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
        [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
        repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))],...
        [min(y) max(y)],'color','r','linestyle','--','linewidth',1);hold on;
    xlim([1 size(A,1)]) 
end


