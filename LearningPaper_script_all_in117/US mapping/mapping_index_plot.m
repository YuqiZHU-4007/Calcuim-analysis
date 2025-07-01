%% overlap of CSmotor 和 CR
figure,
for ii=37:42
    subplot(3,2,ii-trial.test(2)+1)
    a=[0 0.05];t=[1:frame.per_cycle]*fs.ca;
    patch1=patch(([frame.cs_start frame.cs_end frame.cs_end frame.cs_start]-1)*fs.ca,...
        [min(a) min(a) max(a) max(a)],...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    plot(t,squeeze(mean(A_r(:,ii,CSmotor_related_tst_ind),3)),'linewidth',1.5);hold on;
    plot(t,squeeze(mean(A_r(:,ii,CS_related_tst_ind),3)),'linewidth',1.5);hold on;
    bb=re_startpoint(find(re_startpoint(:,1)==ii & re_startpoint(:,2) < frameb.cs_end & re_startpoint(:,2) >= frameb.cs_start-1),2)/fs.behavior;
    plot([bb bb],a,'k--','linewidth',1.5);
    title(['Trial ' num2str(ii)])
end
legend('CS','CS motor related ind in tst','CS related ind in tst')
%% CR in hab
ind=find(CR_ind_up(:,1)==1);
ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);
ind2=intersect(ind2,CS_related_hab_down_ind);
ind3=union(ind,ind2);
ind3=setdiff(1:size(A,2),ind3);
figure,
a=[0 0.03];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
trial_ind=[trial.hab(2):trial.hab(3)];
a=mean(mean(A_r(:,trial_ind,ind3),3),2);plot(a,'color',[0.5 0.5 0.5],'linewidth',2)
b=mean(mean(A_r(:,trial_ind,ind),3),2);plot(b,'color','r','linewidth',2)
c=mean(mean(A_r(:,trial_ind,ind2),3),2);plot(c,'color','b','linewidth',2)
[length(ind3),length(ind2),length(ind)]./size(A,2)

figure,
scatter3(supervoxel(ind3,1),supervoxel(ind3,2),supervoxel(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'b','filled');hold on
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'r','filled');hold on
axis equal;legend('no CR in hab','down CR in hab','up CR in hab')
view(0,90);
%% CR in test
ind=find(CR_ind_up(:,2)==1);
ind=intersect(ind,CS_related_tst_ind);%ind=setdiff(ind,CSmotor_related_tst_ind );
ind2=find(CR_ind_down(:,2)==1);
ind2=intersect(ind2, CS_related_tst_down_ind);
ind3=union(ind,ind2);
ind3=setdiff(1:size(A,2),ind3);
figure,
a=[0 0.03];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
trial_ind=[trial.test(2):trial.test(3)];
a=mean(mean(A_r(:,trial_ind,ind3),3),2);plot(a,'color',[0.5 0.5 0.5],'linewidth',2)
b=mean(mean(A_r(:,trial_ind,ind),3),2);plot(b,'color','r','linewidth',2)
c=mean(mean(A_r(:,trial_ind,ind2),3),2);plot(c,'color','b','linewidth',2)
[length(ind3),length(ind2),length(ind)]./size(A,2)

figure,
scatter3(supervoxel(ind3,1),supervoxel(ind3,2),supervoxel(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'b','filled');hold on
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'r','filled');hold on
axis equal;legend('no CR in tst','down CR in tst','up CR in tst')
view(0,90);
%% test期间 CS引起视觉反应
figure,
a=[0 0.05];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
plot(squeeze(mean(A_r(:,[trial.test(2):trial.test(2)+5],CS_related_tst_ind),3)),'linewidth',2)%trial.test(2):trial.test(3)
figure,
ind=find(sum(UR_ind(:,:),2)==trial.acq_block_num);
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'b','filled');hold on
scatter3(supervoxel(CS_related_tst_ind,1),supervoxel(CS_related_tst_ind,2),supervoxel(CS_related_tst_ind,3),10,'r','filled');
axis equal;legend('UR ind','CS related ind in tst')
%% hab期间 CS引起视觉反应
figure,
a=[0 0.05];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
plot(squeeze(mean(A_r(:,[trial.test(2):trial.test(2)+5],CS_related_hab_ind),3)))%trial.test(2):trial.test(3)
figure,
scatter3(supervoxel(CS_related_tst_ind,1),supervoxel(CS_related_tst_ind,2),supervoxel(CS_related_tst_ind,3),10,'r','filled');hold on
scatter3(supervoxel(CS_related_hab_ind,1),supervoxel(CS_related_hab_ind,2),supervoxel(CS_related_hab_ind,3),10,'b','filled');
axis equal;legend('CS related tst ind','CS related ind in hab')
%% US引起尾动
figure,
ind=find(sum(UR_ind(:,:),2)==trial.acq_block_num);
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'b','filled');hold on
scatter3(supervoxel(USmotor_related_ind,1),supervoxel(USmotor_related_ind,2),supervoxel(USmotor_related_ind,3),10,'r','filled');
axis equal;legend('UR ind','USmotor_related_ind')
%% overlap of CR in test & UR
ind=find(sum(UR_ind(:,:),2)==trial.acq_block_num);ind=setdiff(ind,USmotor_related_ind);
inter=intersect(ind,CS_related_tst_ind);
figure,
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'b','filled');hold on
scatter3(supervoxel(CS_related_tst_ind,1),supervoxel(CS_related_tst_ind,2),supervoxel(CS_related_tst_ind,3),10,'r','filled');
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),10,'k','filled');
axis equal;legend('UR ind non motor','CS related tst ind','overlap')
%% overlap of CRmotor in test & URmotor
ind=USmotor_related_ind;ind2=CSmotor_related_tst_ind;
inter=intersect(ind,ind2);
figure,
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'b','filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'r','filled');
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),10,'k','filled');
axis equal;legend('USmotor related ind','CSmotor related tst ind','overlap')
%%  non-overlap of hab CR and test CR
ind=CS_related_tst_ind;ind2=CS_related_hab_ind;
inter=setdiff(ind,ind2);
figure,
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'b','filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'r','filled');
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),10,'k','filled');
axis equal;legend('CS related tst ind','CS related hab ind','non-overlap')
%% down-regulated CR
ind=find(CR_ind_down(:,2)==1);ind2=find(CR_ind_down(:,1)==1);
inter=setdiff(ind,ind2);
figure,
scatter3(supervoxel(ind,1),supervoxel(ind,2),supervoxel(ind,3),10,'r','filled');hold on
scatter3(supervoxel(ind2,1),supervoxel(ind2,2),supervoxel(ind2,3),10,'b','filled');
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),5,'k','filled');
axis equal;legend('down CR in tst','down CR in hab','non-overlap')
figure,
a=[0 0.05];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
plot(mean(mean(A_r(:,[trial.hab(2):trial.hab(2)+5],ind),3),2))%trial.test(2):trial.test(3)

%% check singal trace 
ind_cluster_regressor_checked_cut=[];ind_cluster_regressor_checked_cut_g=[];
for ii=1:length(ind_cluster_regressor_checked)
    ind_cluster_regressor_checked_cut=[ind_cluster_regressor_checked_cut;ind_cluster_regressor_checked{ii}];
    ind_cluster_regressor_checked_cut_g=[ind_cluster_regressor_checked_cut_g;ii*ones(length(ind_cluster_regressor_checked{ii}),1)];
end

for ii=1:length(ind_cluster_regressor_checked)
iddd=ind_cluster_regressor_checked{ii};clr=GetColormap('hsv_new',length(iddd));
figure,sepplot(1:size(A,1),A(:,iddd),clr);hold on;
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

a=[];
M=A_r(frame_ind,:,:);M=reshape(M,length(frame_ind)*trial.total,[])';
pushbutton_popupplot_Callback(M,ind_cluster_regressor_checked_cut,ind_cluster_regressor_checked_cut_g,...
    GetColormap('hsv_new',length(ind_cluster_regressor_checked)),[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
for ii=1:length(ind_cluster_regressor_checked)
    a(:,ii)=mean(A(:,ind_cluster_regressor_checked{ii}),2);
    length(find(gIX_1toX==ii))
end
figure,sepplot(1:size(A,1),a,GetColormap('hsv_new',length(ind_cluster_regressor_checked)));hold on;
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
clr=GetColormap('hsv_new',length(ind_cluster_regressor_checked));
for ii=1:29
    scatter3(supervoxel(ind_cluster_regressor_checked{ii},1),supervoxel(ind_cluster_regressor_checked{ii},2),supervoxel(ind_cluster_regressor_checked{ii},3),10,clr(ii,:),'filled');hold on;axis equal;
end

%% 统计脑区占比
%hab
ind=find(CR_ind_up(:,1)==1);
ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);
ind2=intersect(ind2,CS_related_hab_down_ind);
ind3=union(ind,ind2);
ind3=setdiff(1:size(A,2),ind3)';
clr=GetColormap(clrmap_name,3);
y=[1*ones(length(ind),1);2*ones(length(ind2),1);3*ones(length(ind3),1)];
x=[ind;ind2;ind3];
load(fullfile(p,'brain arearegion_mask.mat'))
%load('E:\A_Data_lightsheet\Data_huc\20190609\fish3\brain arearegion_mask.mat')
[~,fraction_in_region_in_clust_c,loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],clr,true);
title(['Fish ' num2str(1,'%02d')]);
%test
ind=find(CR_ind_up(:,2)==1);
ind=intersect(ind,CS_related_tst_ind);%ind=setdiff(ind,CSmotor_related_tst_ind );
ind2=find(CR_ind_down(:,2)==1);
ind2=intersect(ind2, CS_related_test_down_ind);
ind3=union(ind,ind2);
ind3=setdiff(1:size(A,2),ind3)';
clr=GetColormap(clrmap_name,3);
y=[1*ones(length(ind),1);2*ones(length(ind2),1);3*ones(length(ind3),1)];
x=[ind;ind2;ind3];
%load(fullfile(p,'brain arearegion_mask.mat'))
[~,fraction_in_region_in_clust_c,loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],clr,true);
title(['Fish ' num2str(1,'%02d')]);
%% mapback
ind=find(CR_ind_up(:,2)==1);
ind=intersect(ind,CS_related_tst_ind);
im=double(env.vol(:,:,1:end-1));rescalegd(im2double(env.vol), [1/10000 1/10000]);%showspv=zeros(size(im,1),size(im,2),3,size(im,3));
showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
for ii=1
    id=ind;
    clr=[1 0 0];kk=1;
    zi=unique(env.supervoxel(id,3));        
    for zz=zi'
      pti =  id(find(env.supervoxel(id,3)==zz));
      slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1],'color',clr(kk,:));
      showspv(:,:,:,zz)=slice;
    end
    kk=kk+1;
end
seqwrite(showspv, fullfile(p,'mapback'));
