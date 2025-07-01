clc;clear all;close all
set(0,'defaultfigurecolor','w')
path=uigetdir('H:\1.Test US\2.Tail free——Data from 117\')
load(fullfile(path,'activities_aft_process.mat'));
load(fullfile(path,'para.mat'));
load(fullfile(path,'env.mat'));
load(fullfile(path,'CR_ind_summary.mat'))
A=activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:);
A_r=reshape(A,frame.per_cycle,trial.total,[]);
B=y_3sd;%B(end+1:frameb.per_cycle*trial.total)=0;
B_r=reshape(B,frameb.per_cycle,trial.total,[]);
supervoxel=env.supervoxel;
for ii=unique(supervoxel(:,3))'
    id=find(supervoxel(:,3)==ii);
    supervoxel(id,3)=(ii-1)*8/0.66+1;
end
supervoxel(:,1)=env.width-supervoxel(:,1);
v=[90,90];colorCS=[0.5 0.5 0.8];
%% task related
% %changes of corr to CS > thres
% thr=0.4;mean(corr_CS_up(:,1))+3*std(corr_CS_up(:,1),[],1)
% figure,hist(corr_CS_up(:,1),[-1:0.1:1]);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
% a=corr_CS_up>thr;
% CS_related_ind=find(sum(a,2)>=1);
% %density figure
% [relative_density,on]=count_spatial_location_density(supervoxel(CS_related_ind,:),supervoxel(:,:), radium);
%   %figure,scatplot(supervoxel(CS_related_ind,1),supervoxel(CS_related_ind,2),'circles', radium,100,[],[],[]);
% %region disribution
%% CS_new_emerged
ind=find(CR_ind_up(:,1)==1);
ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);
ind2=intersect(ind2,CS_related_hab_down_ind);
ind3=union(ind,ind2);
ind3=setdiff(1:size(CR_ind_up(:,1),1),ind3);
ind2=find(CR_ind_up(:,2)==1);
ind2=intersect(ind2,CS_related_tst_ind);
CS_new_emerged=intersect(ind2,ind3);
%% 选取Cb
if exist(fullfile(path,'mask_Hb.mat'))~=0
    load(fullfile(path,'mask_Hb.mat'));
else
    temp=max(env.vol(:,:,1:24),[],3);
    showspv=[];showspv(:,:,1)=temp;showspv(:,:,2)=temp;showspv(:,:,3)=temp;showspv=showspv/255/2;
    showspv = insertShape(showspv, 'filledcircle', [env.supervoxel(CS_new_emerged, 1), env.supervoxel(CS_new_emerged, 2), floor(env.supervoxel(CS_new_emerged, 5)) - 1],'color','r');
    figure,imshow(showspv,[min(showspv(:)) max(showspv(:))]);hold on;
    mask=getmask_imfreehand(temp,1);
    figure,imshow(mask,[min(mask(:)) max(mask(:))]);hold on;
    save(fullfile(path,'mask_Hb.mat'),'mask');
end
B = bwboundaries(mask, 'noholes');plot(B{1}(:,2),B{1}(:,1));
xv=B{1}(:,2);yv=B{1}(:,1);xq=env.supervoxel(CS_new_emerged, 1);yq=env.supervoxel(CS_new_emerged, 2);
[in,on] = inpolygon(xq,yq,xv,yv );
figure,scatter3(supervoxel(CS_new_emerged, 1),supervoxel(CS_new_emerged, 2),supervoxel(CS_new_emerged, 3));hold on;
scatter3(supervoxel(CS_new_emerged(in), 1),supervoxel(CS_new_emerged(in), 2),supervoxel(CS_new_emerged(in), 3),'filled');
axis equal;view(v);
Cb_ind_new_emerged=CS_new_emerged(in);
length(CS_new_emerged)/size(env.supervoxel,1)
%% Cb内的CR corr
figure('position',[162,558,1645,420]),t=tiledlayout(4,1);t.TileSpacing = 'compact';t.Padding = 'compact';
aa=A_r(:,:,Cb_ind_new_emerged);aa=reshape(aa,frame.per_cycle*trial.total,[]);
%aa=smoothdata(aa,'movmedian',3);
ax1=nexttile(t,[1,1]);plot([1:size(A,1)].*fs.ca,mean(aa,2),'linewidth',2,'color','k');hold on;
y=[-0.01 0.1];
patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']'.*fs.ca,...
    repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
if size(startpoint,2)==1
    line([startpoint,startpoint]./fs.behavior,[min(y) max(y)],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',0.5);hold on;
else
    line([startpoint;startpoint]./fs.behavior,[min(y) max(y)],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',0.5);hold on;
end
line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))].*fs.ca,...
    [min(y) max(y)],'color','r','linestyle','--','linewidth',2);hold on;
xlim([1 size(A,1)].*fs.ca);ylim(y);set(gca,'fontsize',16);box off;xticklabels(ax1,{});

frame_ind=frame.cs_start:frame.us_start-1;trial_ind=1:trial.total;
a=A_r(frame_ind,trial_ind,Cb_ind_new_emerged);a=reshape(a,length(frame_ind)*length(trial_ind),[]);
nexttile(t,[3,1]);imagesc(a',y);colormap('hot');hold on;
%line(repmat([frame.cs_start:length(frame_ind):length(frame_ind)*length(trial_ind)],2,1),repmat([0 length(Cb_ind_new_emerged)],length(trial_ind),1)','color','w','linestyle','--','linewidth',1);

%correlation map
figure('position',[389,644,851,334]),t=tiledlayout(1,2);t.TileSpacing = 'compact';t.Padding = 'compact';
n=2;
Cor=corr(a);
[idx_functional_CR,~,~,~] = kmeans(Cor',n,'Distance','correlation','Replicates',10);
[idxAs_CR,idxAa_ind_CR]=sort(idx_functional_CR);
nexttile(t,[1,1]);imagesc(Cor(idxAa_ind_CR,idxAa_ind_CR),[0 1]);colormap('jet');colorbar;%axis equal
clr=[1.00,0.41,0.16;0.39,0.83,0.07];GetColormap('hsv_new',n);
[h]=pushbutton_popupplot_Callback(a',idxAa_ind_CR,idxAs_CR,clr,[],fs.ca,ones(1,length(frame_ind)*length(trial_ind)),ones(1,length(frame_ind)*length(trial_ind)),1,'all',0,0,0,0);
nexttile(t,[1,1]);
scatter3(supervoxel(CS_new_emerged, 1),supervoxel(CS_new_emerged, 2),supervoxel(CS_new_emerged, 3),10,[0.5 0.5 0.5],'filled');hold on;
for ii=1:n
    ind=Cb_ind_new_emerged(idxAa_ind_CR(find(idxAs_CR==ii) ));
    scatter3(supervoxel(ind, 1),supervoxel(ind, 2),supervoxel(ind,3),12,clr(ii,:),'filled');hold on;grid off;
    axis off;
end
legend;
axis equal;view(v)
%% CR profile & AUC of each cluster
% 按CS onset平均;Go trial & nogo trial CR
AUC_go_clust=[];AUC_nogo_clust=[];AUC_us_clust=[];nogo_clust=[];go_clust=[];p_clust=[];
%y=[-0.01 0.02];
for hh=1:n
    ind=Cb_ind_new_emerged(idxAa_ind_CR(find(idxAs_CR==hh) ));
    figure('position',[28,330,1780,319]),t=tiledlayout(2,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
    for ii=1:trial.acq_block_num+2
        switch ii
            case 1
                trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';b=frame.cs_end-3;bb=frameb.cs_end-10;
            case trial.acq_block_num+2
                trial_ind=[trial.test(2):trial.test(3)];leg='Post';b=frame.cs_end-3;bb=frameb.cs_end-10;
            otherwise
                trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
                b=frame.us_start-3;bb=frameb.us_start-10;
        end
        go=[];kk=1;
        for jj=1:length(trial_ind)
            a=re_startpoint(find(re_startpoint(:,1)==trial_ind(jj) ...
                & re_startpoint(:,2)>=frameb.cs_start-10 & re_startpoint(:,2)<bb),1);
            if ~isempty(a)
                go(kk)=unique(a);kk=kk+1;
            end
        end
        nogo=setdiff(trial_ind,go);
        y1=squeeze(mean(A_r(:,go,ind),2));
        y2=squeeze(mean(A_r(:,nogo,ind),2));
        AUC_go_clust{hh}(:,ii)= sum(y1(frame.cs_start:b,:));
        AUC_nogo_clust{hh}(:,ii)= sum(y2(frame.cs_start:b,:));
        x=(1:frame.per_cycle)*fs.ca;
        nexttile(t,[1 1]);
        patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start]'.*fs.ca,...
            [min(y) min(y) max(y) max(y)],...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        if ii>1 && ii< trial.acq_block_num+2
            line([frame.us_start,frame.us_start].*fs.ca,[min(y) max(y)],'color','k','linestyle','--','linewidth',2);hold on;
        end
        yy=mean(y1,2);%yy(find(yy>0.04))=mean(yy);
        if ~isempty(y1) shadedErrorBar(x,yy,[std(y1,[],2) , std(y1,[],2)],'lineProps','b');hold on;end
        if ~isempty(y1) shadedErrorBar(x,mean(y2,2),[std(y2,[],2) , std(y2,[],2)],'lineProps','r');hold on;end
        set(gca,'fontsize',16);xlim([1 frame.per_cycle]*fs.ca);ylim(y);box on;
        %if ii==trial.acq_block_num+2 legend('CS','US','GO','No GO');end
    end
    % 按行为onset平均
    %figure('position',[28,439,1780,210]),t=tiledlayout(1,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
    for ii=1:trial.acq_block_num+2
        switch ii
            case 1
                trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';b=frame.cs_end-1;bb=frameb.cs_end-1;
            case trial.acq_block_num+2
                trial_ind=[trial.test(2):trial.test(3)];leg='Post';b=frame.cs_end-1;bb=frameb.cs_end-1;
            otherwise
                trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
                b=frame.us_start-1;bb=frameb.us_start-1;
        end
        kk=1;aa=[];
        for jj=1:length(trial_ind)
            a=find(re_startpoint(:,1)==trial_ind(jj) ...
                & re_startpoint(:,2)>=frameb.cs_start-10 & re_startpoint(:,2)<bb);
            if ~isempty(a)
                a=a(1);
                b=round(startpoint(a)/fs.behavior/fs.ca)-ceil(frame.per_cycle/2):min(round(startpoint(a)/fs.behavior/fs.ca)+ceil(frame.per_cycle/2)-1,frame.per_cycle*trial.total);
                aa(1:length(b),kk,:)=A(b, ind);
                kk=kk+1;
            end
        end
        aa=squeeze(mean(aa,2));
        x=(1:ceil(frame.per_cycle/2)*2)*fs.ca;
        nexttile(t,[1 1]);
        line([ceil(frame.per_cycle/2),ceil(frame.per_cycle/2)].*fs.ca,[min(y) max(y)],'color','k','linestyle','--','linewidth',2);hold on;
        if ~isempty(aa) shadedErrorBar(x,mean(aa,2),[std(aa,[],2) , std(aa,[],2)],'lineProps','r');hold on;end
        set(gca,'fontsize',16);xlim([1 frame.per_cycle]*fs.ca);ylim(y);box on;
        %if ii==trial.acq_block_num+2 legend('CS','US','GO','No GO');end
    end
    for ii=1:trial.acq_block_num
        trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-1)*trial.acq_block_trial;
        y1=squeeze(mean(A_r(:,trial_ind,ind),2));
        AUC_us_clust{hh}(:,ii)= sum(y1(frame.us_start:frame.us_start+1,:));
    end
    figure,plot(repmat([1:trial.acq_block_num],length(ind),1),AUC_us_clust{hh},'marker','o');
    go_clust{hh}(:,1)=mean(AUC_go_clust{hh});
    go_clust{hh}(:,2)=std(AUC_go_clust{hh},[],1);
    go_clust{hh}(:,3)=repmat(length(AUC_go_clust{hh}),size(AUC_go_clust{hh},2),1); 
    nogo_clust{hh}(:,1)=mean(AUC_nogo_clust{hh});
    nogo_clust{hh}(:,2)=std(AUC_nogo_clust{hh},[],1);
    nogo_clust{hh}(:,3)=repmat(length(AUC_nogo_clust{hh}),size(AUC_nogo_clust{hh},2),1);
    for ii=1:size(AUC_nogo_clust{hh},2)
        if ~isnan(AUC_nogo_clust{hh}(1,ii))
            p_clust(ii,hh)=signrank(AUC_nogo_clust{hh}(:,1),AUC_nogo_clust{hh}(:,ii),'tail','left');
        end
    end
    find(p_clust(:,hh)<0.05)
end
%% CR profile & AUC
% 按CS onset平均;不区分Go trial & nogo trial CR
AUC=[];
figure('position',[28,439,1780,210]),t=tiledlayout(1,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';b=frame.cs_end-1;bb=frameb.cs_end-10;
        case trial.acq_block_num+2
            trial_ind=[trial.test(2):trial.test(3)];leg='Post';b=frame.cs_end-1;bb=frameb.cs_end-10;
        otherwise
            trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
            b=frame.us_start-1;bb=frameb.us_start-10;
    end
    y1=squeeze(mean(A_r(:,trial_ind,Cb_ind_new_emerged),2));
    AUC(:,ii)= sum(y1(frame.cs_start:b,:));
    x=(1:frame.per_cycle)*fs.ca;
    nexttile(t,[1 1]);
    patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start]'.*fs.ca,...
        [min(y) min(y) max(y) max(y)],...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    if ii>1 && ii< trial.acq_block_num+2
        line([frame.us_start,frame.us_start].*fs.ca,[min(y) max(y)],'color','k','linestyle','--','linewidth',2);hold on;
    end
    if ~isempty(y1) shadedErrorBar(x,mean(y1,2),[std(y1,[],2) , std(y1,[],2)],'lineProps','r');hold on;end
    set(gca,'fontsize',16);xlim([1 frame.per_cycle]*fs.ca);ylim(y);box on;
    %if ii==trial.acq_block_num+2 legend('CS','US','GO','No GO');end
end
all=[];all(:,1)=mean(AUC);
all(:,2)=std(AUC,[],1);
all(:,3)=repmat(length(AUC),size(AUC,2),1);
p=[];
for ii=1:size(AUC,2)
    if ~isnan(AUC(1,ii))
    p(ii)=signrank(AUC(:,1),AUC(:,ii),'tail','left');
    end
end
find(p<0.05)
% 按CS onset平均;Go trial & nogo trial CR
AUC_go=[];AUC_nogo=[];AUC_us=[];
figure('position',[28,439,1780,210]),t=tiledlayout(1,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';b=frame.cs_end-1;bb=frameb.cs_end-10;
        case trial.acq_block_num+2
            trial_ind=[trial.test(2):trial.test(3)];leg='Post';b=frame.cs_end-1;bb=frameb.cs_end-10;
        otherwise
            trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
            b=frame.us_start-1;bb=frameb.us_start-10;
    end
    go=[];kk=1;
    for jj=1:length(trial_ind)
        a=re_startpoint(find(re_startpoint(:,1)==trial_ind(jj) ...
            & re_startpoint(:,2)>=frameb.cs_start-10 & re_startpoint(:,2)<bb),1);
        if ~isempty(a)
            go(kk)=unique(a);kk=kk+1;
        end
    end
    nogo=setdiff(trial_ind,go);
    y1=squeeze(mean(A_r(:,go,Cb_ind_new_emerged),2));
    y2=squeeze(mean(A_r(:,nogo,Cb_ind_new_emerged),2));
    AUC_go(:,ii)= sum(y1(frame.cs_start:b,:));
    AUC_nogo(:,ii)= sum(y2(frame.cs_start:b,:));
    x=(1:frame.per_cycle)*fs.ca;
    nexttile(t,[1 1]);
    patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start]'.*fs.ca,...
        [min(y) min(y) max(y) max(y)],...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    if ii>1 && ii< trial.acq_block_num+2
        line([frame.us_start,frame.us_start].*fs.ca,[min(y) max(y)],'color','k','linestyle','--','linewidth',2);hold on;
    end
    if ~isempty(y1) shadedErrorBar(x,mean(y1,2),[std(y1,[],2) , std(y1,[],2)],'lineProps','b');hold on;end
    if ~isempty(y1) shadedErrorBar(x,mean(y2,2),[std(y2,[],2) , std(y2,[],2)],'lineProps','r');hold on;end
    set(gca,'fontsize',16);xlim([1 frame.per_cycle]*fs.ca);ylim(y);box on;
    %if ii==trial.acq_block_num+2 legend('CS','US','GO','No GO');end
end
% 按行为onset平均
figure('position',[28,439,1780,210]),t=tiledlayout(1,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';b=frame.cs_end-1;bb=frameb.cs_end-1;
        case trial.acq_block_num+2
            trial_ind=[trial.test(2):trial.test(3)];leg='Post';b=frame.cs_end-1;bb=frameb.cs_end-1;
        otherwise
            trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
            b=frame.us_start-1;bb=frameb.us_start-1;
    end
   kk=1;aa=[];
    for jj=1:length(trial_ind)
        a=find(re_startpoint(:,1)==trial_ind(jj) ...
            & re_startpoint(:,2)>=frameb.cs_start-10 & re_startpoint(:,2)<bb);
        if ~isempty(a)        
            a=a(1);
            b=round(startpoint(a)/fs.behavior/fs.ca)-ceil(frame.per_cycle/2):min(round(startpoint(a)/fs.behavior/fs.ca)+ceil(frame.per_cycle/2)-1,frame.per_cycle*trial.total);
            aa(1:length(b),kk,:)=A(b, Cb_ind_new_emerged);
            kk=kk+1;
        end
    end
    aa=squeeze(mean(aa,2));
    x=(1:ceil(frame.per_cycle/2)*2)*fs.ca;
    nexttile(t,[1 1]);
    line([ceil(frame.per_cycle/2),ceil(frame.per_cycle/2)].*fs.ca,[min(y) max(y)],'color','k','linestyle','--','linewidth',2);hold on;
    if ~isempty(aa) shadedErrorBar(x,mean(aa,2),[std(aa,[],2) , std(aa,[],2)],'lineProps','r');hold on;end
    set(gca,'fontsize',16);xlim([1 frame.per_cycle]*fs.ca);ylim(y);box on;
    %if ii==trial.acq_block_num+2 legend('CS','US','GO','No GO');end
end
for ii=1:trial.acq_block_num
    trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-1)*trial.acq_block_trial;
    y1=squeeze(mean(A_r(:,trial_ind,Cb_ind_new_emerged),2));
    AUC_us(:,ii)= max(y1(frame.us_start:frame.us_start+2,:))-mean(y1(frame.cs_start:frame.us_start-1,:));
end
figure,plot(repmat([1:trial.acq_block_num],length(Cb_ind_new_emerged),1),AUC_us,'marker','o');
% auto k-means
%统计每个阶段有CR的比例
response=[];unresponse=[];
figure('position',[389,644,851,334]),t=tiledlayout(1,trial.acq_block_num);t.TileSpacing = 'compact';t.Padding = 'compact';
for ii=1:trial.acq_block_num
    response(:,ii)=length(find( CR_ind_up_acq(Cb_ind_new_emerged,ii)==1));
    unresponse(:,ii)=length(find( CR_ind_up_acq(Cb_ind_new_emerged,ii)~=1));
    nexttile(t,[1,1]);pie([response(:,ii) unresponse(:,ii)])
end
legend('CS Responsive','CS UnResponsive');
go=[];go(:,1)=mean(AUC_go);
go(:,2)=std(AUC_go,[],1);
go(:,3)=repmat(length(AUC_go),size(AUC_go,2),1);

nogo=[];nogo(:,1)=mean(AUC_nogo);
nogo(:,2)=std(AUC_nogo,[],1);
nogo(:,3)=repmat(length(AUC_nogo),size(AUC_nogo,2),1);
p=[];
for ii=1:size(AUC_nogo,2)
    if ~isnan(AUC_nogo(1,ii))
    p(ii)=signrank(AUC_nogo(:,1),AUC_nogo(:,ii),'tail','left');
    end
end
find(p<0.05)
response_al=[length(find( CR_ind_up(Cb_ind_new_emerged,1)==1)),response,length(find( CR_ind_up(Cb_ind_new_emerged,2)==1))]
a=(response_al/length(Cb_ind_new_emerged))';
IPower=all(:,1).*(response_al/length(Cb_ind_new_emerged))';
save(fullfile(path,'Hb_summary.mat'),'AUC','AUC_go','AUC_nogo','AUC_us',...
    'idx_functional_CR','AUC_go_clust','AUC_nogo_clust','AUC_us_clust','response_al','IPower');
%% 自发活动
if 0
if exist(fullfile(path,'act_spon.mat'))
    load(fullfile(path,'act_spon.mat'));
    %dfdf
    dfdf_spon=zeros(size(activities_spon));
    for ii=1:size(activities_spon,3)
        winsize=41;
        activities=squeeze(activities_spon(:,:,ii));
        activities=smoothdata(activities,2,'movmedian',3);
        activities=normalize(activities,2,'zscore');
        [a,b]=compute_dfdf(activities,winsize);
        dfdf_spon(:,:,ii)=a;
        % [a,b]=compute_dfdf(mean(activities(CS_new_emerged,:),1),winsize);
        % figure,plot(mean(activities(CS_new_emerged,:),1));hold on;plot(b)
        %figure,plot(activities(CS_new_emerged(100),:));hold on;plot(b(CS_new_emerged(100),:))
    end
    %% Cb和其他的Cor
    %自发
    %correlation map（所有emerged 和 Cb内部的）
    for zz=1:2
        switch zz
            case 1
                ind=CS_new_emerged;
            case 2
                ind=Cb_ind_new_emerged;
        end
        n=3;figure('position',[299,522,1624,456]);
        for ii=size(activities_spon,3):-1:1
            a=squeeze(dfdf_spon(ind,:,ii));
            Cor=corr(a');
            if ii==size(activities_spon,3)
                %figure,imagesc(Cor,[0 max(Cor(:))]);colormap('hot');colorbar;
                [idx_functional,~,~,~] = kmeans(Cor',n,'Distance','correlation','Replicates',10);
                [idxAs,idxAa_ind]=sort(idx_functional);
            end
            %figure;hist(idx_functional,[1:n]);xlim([1 n]);
            subplot(2,5,ii);imagesc(Cor(idxAa_ind,idxAa_ind),[-0.2 0.2]);colormap('jet');colorbar;%axis equal
            title(num2str(ii))
        end
        clr=GetColormap('hsv_new',n);figure;
        for ii=1:n
            indd=ind(idxAa_ind(find(idxAs==ii) ));
            scatter3(supervoxel(indd, 1),supervoxel(indd, 2),supervoxel(indd,3),12,clr(ii,:),'filled');hold on;
        end
        axis equal;view(v)
    end
    %按CR cor cluster 排序自发活动cor
    figure('position',[299,522,1624,456]);
    [~,idxAa_ind]=sort(idx_functional_CR);
    for ii=size(activities_spon,3):-1:1
        a=squeeze(dfdf_spon(Cb_ind_new_emerged,:,ii));
        Cor=corr(a');
        subplot(2,5,ii);imagesc(Cor(idxAa_ind,idxAa_ind),[0 0.2]);colormap('jet');colorbar;
        title(num2str(ii))
    end
end
end
%CS
%% Cb的dpca

%% save
if ~isempty(Cb_ind_new_emerged)
radium=floor(15/0.66);
[Spatial_density,On]=count_spatial_location_density(supervoxel(Cb_ind_new_emerged,:),supervoxel(:,:), radium);
length(Cb_ind_new_emerged)
mean(Spatial_density)
else
    Spatial_density=[];
    idxAa_ind_CR=[];idxAs_CR=[];
end
save(fullfile(path,'mask_Hb.mat'),'mask','Cb_ind_new_emerged','Spatial_density','AUC_go','AUC_nogo','AUC_us',...
'idxAa_ind_CR','idxAs_CR','AUC_go_clust','AUC_nogo_clust','AUC_us_clust');

