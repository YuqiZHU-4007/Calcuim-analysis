%% 0 to 1 CR 
ind=find(CR_ind_up(:,1)==1);ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);ind2=intersect(ind2, CS_related_hab_down_ind);
ind3=union(ind,ind2);ind3=setdiff(1:size(A,2),ind3);
ind=find(CR_ind_up(:,2)==1);
ind=intersect(ind,CS_related_tst_ind);%ind=setdiff(ind,CSmotor_related_tst_ind );
inter=intersect(ind,ind3);
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
numK=[];[gIXx,~,numK,~]=find_best_k_in_range(M,3:50);
cIXx=1:length(inter);cIX_reg = (1:size(M,1))';
masterthres=0.6;
para=struct;para=setfield(para,'k1',numK);
para=setfield(para,'merge',masterthres);
para=setfield(para,'cap',masterthres);
para=setfield(para,'reg1',masterthres);
para=setfield(para,'reg2',masterthres);
para=setfield(para,'minSize',5);
for ii=1
[cIX_0to1,gIX_0to1] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX_0to1))
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_0to1));
vi='on';
for iii=1
    cIX_iii=cIX_0to1;gIX_iii=gIX_0to1;clrmap_iii=clrmap;
    %     [h,gIX2, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M)
    %     length(unique(gIX2))
    %2
    B=M;
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_0to1))
        a(:,ii)=mean(A(:,inter(cIX_0to1(find(gIX_0to1==ii)))),2);
        length(find(gIX_0to1==ii))
    end
    figure,sepplot(1:size(A,1),a,GetColormap('hsv_new',max(gIX_0to1)));hold on;
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
    indd_0to1={};aa=[];indd_0to1_regressor={};
    for ii=1:length(unique(gIX_0to1))
        M_all=A_r(frame_ind,:,:);M_all=reshape(M_all,length(frame_ind)*trial.total,[]);
        [stimcorr,~] = MotorSourceCorrelation(M_all',mean(M(cIX_0to1(find(gIX_0to1==ii)),:),1),[]);
        n=3;
        thr=max(mean(stimcorr)+3.5*std(stimcorr),0.5)
        %figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
        indd_0to1{ii}=intersect(find(stimcorr>thr), ind3);
        indd_0to1_regressor{ii}=inter(cIX_0to1(find(gIX_0to1==ii)));
        aa(:,ii)=mean(A(:,indd_0to1{ii}),2);
        %figure,plot(a(:,ii));hold on;plot(mean(A(:,indd{ii}),2),'r');legend('regressor','stim output')
        subplot(4,4,ii),
        scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),2,[0.5 0.5 0.5],'filled');hold on
        scatter3(supervoxel(indd_0to1{ii},1),supervoxel(indd_0to1{ii},2),supervoxel(indd_0to1{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
        title(num2str(ii));%legend;
        view(-90,90)
    end
    figure,
    for ii=1:length(unique(gIX_0to1))
        scatter3(supervoxel(indd_0to1_regressor{ii},1),supervoxel(indd_0to1_regressor{ii},2),supervoxel(indd_0to1_regressor{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    figure,sepplot([1:size(A,1)].*fs.ca,aa,GetColormap('hsv_new',max(gIX_0to1)));hold on;
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
%% cluster间的correlation
a=[];kk=1;
for ii=1:length(unique(gIX_0to1))
    a(:,kk:kk+length(find(gIX_0to1==ii))-1)=A(:,inter(cIX_0to1(find(gIX_0to1==ii))));
    kk=kk+length(find(gIX_0to1==ii));
end
c=corr(a);
figure,imagesc(c,[0 0.5]);colorbar;

%% go-trial 和 nogo trial
go_ind=[37:42];[37 38 40 41 42];nogo_ind=[39 43 44 45 46];
clr=GetColormap('hsv_new',max(gIX_0to1));
%按CS对齐
figure,subplot(3,1,1)
for ii=1:length(unique(gIX_0to1))
    a=mean(squeeze(mean(A_r(:,go_ind,indd_0to1{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('Go trial');%legend;
subplot(3,1,2)
for ii=1:length(unique(gIX_0to1))
    a=mean(squeeze(mean(A_r(:,nogo_ind,indd_0to1{ii}),3)),2);
    plot(a,'color',clr(ii,:),'linewidth',2);hold on;
end
title('NoGo trial');%legend;
%按运动对齐
subplot(3,1,3)
for ii=1:length(unique(gIX_0to1))
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
%% 统计脑区占比
x=[];y=[];
for ii=1:length(unique(gIX_0to1))
    x=[x,indd_0to1{ii}];
    y=[y,ii*ones(length(indd_0to1{ii}),1)'];
end
y=[y,(max(y)+1)*ones(length(setdiff(1:size(A,2),x)),1)']';
x=[x setdiff(1:size(A,2),x)];
clr=GetColormap('hsv_new',length(unique(gIX_0to1))+1);
load('E:\A_Data_lightsheet\Data_huc\20190609\fish3\brain arearegion_mask.mat')
[~,fraction_in_region_in_clust_c,loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],clr,true);
title(['Fish ' num2str(1,'%02d')]);
