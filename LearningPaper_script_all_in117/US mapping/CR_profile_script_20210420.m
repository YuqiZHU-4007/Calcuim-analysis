%% 计算CR 面积
set(0,'defaultfigurecolor','w');
iddd=[indd_1toX,indd_0to1];iddd{10}=[iddd{9} iddd{10}];
index_cluster=[2 5 10 11 7 8];clr=GetColormap('hsv_new',length(index_cluster));
win=frame.cs_start:frame.us_start-1;
aa=[];indd_all=[];indd_all_g=[];
for ii=1:length(index_cluster)
    a=squeeze(sum(A_r(win,:,iddd{index_cluster(ii)}),1));
    errBar=std(a,[],1);
    y=mean(a,1);x=1:trial.total;
    figure('position',[284,674,530,137]), shadedErrorBar(x,a',{@mean,@std},'lineprops',{clr(ii,:)}); hold on;
    scatter(x,mean(a,2),'MarkerFaceColor',clr(ii,:))  
    xlim([1 trial.total]);
    set(gca,'fontsize',16);box on;
    aa(:,ii)=mean(A(:,iddd{index_cluster(ii)}),2);
    indd_all=cat(2,indd_all,iddd{index_cluster(ii)});
    indd_all_g=[indd_all_g,ii*ones(1,length(iddd{index_cluster(ii)}))];
end
bb={}; 
figure('position',[480,405,230*4,273*2]),
for ii=1:length(index_cluster)
    a=squeeze(sum(A_r(win,:,iddd{index_cluster(ii)}),1));
    x=mean(a(1:trial.hab(1),:),1);y=mean(a(trial.test(2):trial.test(3),:),1);
    p=signrank(x,y,'tail','left')
    bb{ii}=[x;y];
    c=categorical({'Pre Cond.','Post Cond.'});
    subplot(2,4,ii),boxplot([x;y]','Labels',{'Pre Cond.','Post Cond.'})
%     figure, bar(c,[mean(x);mean(y)]'); hold on;
%     scatter(repmat(c(1),1,length(x)),x)
    set(gca,'fontsize',14);box on;
end

figure,sepplot(1:size(A,1),aa,clr);hold on;
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
M=A_r(frame_ind,:,:);M=reshape(M,length(frame_ind)*trial.total,[])';M_norm=normalize(M,2,'zscore');M=M_norm;
[h]=pushbutton_popupplot_Callback(M,indd_all,indd_all_g,clr,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);

for ii=1:length(index_cluster)
       figure,
   scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),2,[0.5 0.5 0.5],'filled');hold on

   scatter3(supervoxel(iddd{index_cluster(ii)},1),supervoxel(iddd{index_cluster(ii)},2),supervoxel(iddd{index_cluster(ii)},3),10,clr(ii,:,:),'filled');hold on;axis equal;
end
legend;


im=double(env.vol(:,:,1:end-1));rescalegd(im2double(env.vol), [1/10000 1/10000]);%showspv=zeros(size(im,1),size(im,2),3,size(im,3));
showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
kk=1;
for ii=1:length(index_cluster)
    id=iddd{index_cluster(ii)};
    zi=unique(env.supervoxel(id,3));        
    for zz=zi'
      pti =  id(find(env.supervoxel(id,3)==zz));
      slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1],'color',clr(ii,:));
      showspv(:,:,:,zz)=slice;
    end
end
seqwrite(showspv, 'F:\DUlab\Progress report\202103\all');

%% correlation with motor regressor
m=[];
for ii=1:length(index_cluster)
    m(ii,1)=mean(mtRescorr_hab(iddd{index_cluster(ii)}));
    m(ii,2)=std(mtRescorr_hab(iddd{index_cluster(ii)}));
    m(ii,3)=length(mtRescorr_hab(iddd{index_cluster(ii)}));
end

m=nan(450,7);
for ii=1:length(index_cluster)
    m(1:length(iddd{index_cluster(ii)}),ii)=mtRescorr_all(iddd{index_cluster(ii)});
end
x=mtRescorr_all(iddd{index_cluster(1)});y=mtRescorr_all(iddd{index_cluster(6)});
p=ranksum(x,y,'tail','left')


%% 不同类的脑区占比
y=[indd_all_g,(max(indd_all_g)+1)*ones(length(setdiff(1:size(A,2),indd_all)),1)']';
x=[indd_all setdiff(1:size(A,2),indd_all)];
clr=GetColormap('hsv_new',length(index_cluster)+1);
load('E:\A_Data_lightsheet\Data_huc\20190605\fish3\brain arearegion_mask.mat')
[~,fraction_in_region_in_clust_c,loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,y,x,...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],...
    [env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],clr,true);
title(['Fish ' num2str(1,'%02d')]);
