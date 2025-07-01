%CR profile(isevent & getpara_event) & classification
%zyq
%behaviorpath='E:\A_Data_lightsheet\Data_huc\20180509_huc\fish2\20180509_fish2_huc';
[actname,actpath]=uigetfile('','act_aftprocess.mat');
[envname,envpath]=uigetfile(actpath,'env.mat');
[bahname,behpath]=uigetfile(actpath,'behav');
[paraname,parapath]=uigetfile(actpath,'para.mat');

load([actpath actname]);
load([envpath envname]);
load([parapath paraname]);

%[activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf_huc(activities');

%[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([behpath bahname]);

trial_ind=trial.hab(2):min(trial.hab(2)+6-1,trial.hab(3)); %%trial.test(2):trial.test(2)+2
ref_win=ceil(0/fs.ca+1):frame.cs_start-1;%取时间段的baseline作为参考判断event
area_win_hab=frame.cs_start:frame.cs_end-1;

act=activities_preCS_dfdf_aftcorrect';
act_CS_hab= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,[],size(act,1));
%act_CS_hab_mean=zeros(frame.per_cycle,size(act,1));
act_CS_hab_mean=mean(act_CS_hab,2);

for ii=1:size(activities_preCS_dfdf_aftcorrect,2)
    a=act_CS_hab_mean(:,:,ii);
    event(1,ii)=isevent_20190508(a,ref_win,area_win_hab,trial);
end

CR_loc=find(event==1); length(CR_loc)
%showspv=mapback_to_vol_zyq_20180520(env,opt,CR_loc,'circle');
showspv= mapback(ones(size(CR_loc)), env.supervoxel,[env.height env.width env.depth],CR_loc);
show_spv_GUI(showspv);

%classication
act_CS_hab_mean_isevent=act_CS_hab_mean(:,CR_loc);
figure,
line([frame.cs_start frame.cs_start],[0 0.015],'color','b','linewidth',1.5,'linestyle','--');hold on;
line([frame.cs_end frame.cs_end],[0 0.015],'color','b','linewidth',1.5,'linestyle','--');hold on;
plot(act_CS_hab_mean_isevent,'linewidth',2,'color',[0.5 0.5 0.5]);hold on
plot(mean(act_CS_hab_mean_isevent,2),'linewidth',3,'color',[1 0 0])

figure,
line([frame.cs_start frame.cs_start],[0 0.015],'color','b','linewidth',1.5,'linestyle','--');hold on;
line([frame.cs_end frame.cs_end],[0 0.015],'color','b','linewidth',1.5,'linestyle','--');hold on;
plot(mean(act_CS_hab_mean_isevent,2),'linewidth',3,'color',[1 0 0]);hold on
para_mean=getpara_event_20190521(mean(act_CS_hab_mean_isevent,2),ref_win,frame.cs_start:frame.per_cycle);
figure,txt_dis=string;
for ii=1:length(name)
    a = getfield(para_mean,name{ii});
    txt_dis=[txt_dis;name{ii},': ',num2str(a)];
end
txt = uicontrol('Style','text',...
    'Position',[10+5+(i-1)*100,50,100,200],...
    'fontsize',10,...
    'BackgroundColor','w',...
    'horizontalalignment','left',...
    'string',txt_dis);

para=getpara_event_20190521(act_CS_hab_mean_isevent,ref_win,frame.cs_start:frame.per_cycle);

name=fieldnames(para);
for ii=1:length(name)
    a = getfield(para,name{ii});
    %a=mynorm(a);
    a=normalize(a,'zscore');
    norm_para(:,ii)=a;
end

norm_para_pca=pca(norm_para');
eucD = pdist(norm_para_pca(:,1:3),'Euclidean');
eucD2 = pdist(norm_para_pca(:,1:3),'cosine');


figure,
subplot(1,2,1),
clustTree = linkage(eucD,'average');
cophenet(clustTree,eucD)
[h,nodes] = dendrogram(clustTree,200);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
subplot(1,2,2),
clustTree = linkage(eucD2,'average');
cophenet(clustTree,eucD2)
[h,nodes] = dendrogram(clustTree,200);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
%%%%%%%%%%%%%%%%%%%k-means
classnum=3;
[cidx3,cmeans3,sumd3] = kmeans(norm_para_pca(:,1:3),classnum,'replicates',100,'display','final','Distance','sqEuclidean');%'cosine','sqEuclidean'
figure,[silh3,h] = silhouette(norm_para_pca(:,1:3),cidx3,'sqEuclidean');
B=act_CS_hab_mean_isevent;
x=1:size(B,1);
col=3;
y=[0 0.02];
figure,    
for i = 1:length(unique(cidx3))
    subplot(ceil(classnum/col),col,i)
    clust = find(cidx3==i);
    %     for jj=1:length(unique(ind))
    %         fishind=find(ind==jj);
    %         inter_ind=intersect(clust,fishind);
    %         text(5,0.1+0.05*jj,[num2str(jj) ': ' num2str(length(inter_ind))],'fontsize',25);hold on;
    %     end
    n(:,i)=length(clust);
    m(:,i)=mean(B(:,clust),2);
    sd(:,i)=std(B(:,clust),[],2);
    plot(x,B(:,clust),'color',[0.5 0.5 0.5],...
        'LineWidth',1);
    hold on
    plot(x,mean(B(:,clust),2),'-','color',[1 0 0],...
        'LineWidth',2);
    hold on
    set(gca,'FontSize',20);box off;xlim([x(1) x(end)]);%ylim([0 1.5]);
    %text(ceil(x/2),0.4,num2str(length(clust)),'fontsize',25);hold on;
    title(['class ' num2str(i) '_' num2str(length(clust))],'fontsize',15,'interpreter','none')
    ylabel('△F/F0');xlabel('frame');
    ylim(y);
end
hold off
y=[0 0.02];
figure, 
para_mean=struct;
for i = 1:length(unique(cidx3))
    clust = find(cidx3==i);
    plot(x,mean(B(:,clust),2),'-',...
        'LineWidth',2);hold on
    para_mean(i)=getpara_event_20190521(mean(B(:,clust),2),ref_win,frame.cs_start:frame.per_cycle);
    set(gca,'FontSize',20);box off;xlim([x(1) x(end)]);%ylim([0 1.5]);
    %text(ceil(x/2),0.4,num2str(length(clust)),'fontsize',25);hold on;
    ylabel('△F/F0');xlabel('frame');
    ylim(y);
end
line([frame.cs_start frame.cs_start],y,'color','b','linewidth',1.5,'linestyle','--');hold on;
line([frame.cs_end frame.cs_end],y,'color','b','linewidth',1.5,'linestyle','--');hold on;
legend({num2str([1:classnum]')})
hold off
y=[0 0.02];
figure,
showspv_mean_class={};
for i = 1:length(unique(cidx3))
    s=subplot(ceil(classnum/col),col,i);
    clust = find(cidx3==i);
    %showspv{i}= mapback(ones(size(clust)), env.supervoxel,[env.height env.width env.depth],CR_loc(clust));
    showspv_mean_class{i}=mapback_to_vol_zyq_20180520(env,opt,CR_loc(clust),'circle');
    plot(x,mean(B(:,clust),2),'-','color',[1 0 0],...
        'LineWidth',2);hold on
    set(gca,'FontSize',20,'visible','off');box off;xlim([x(1) x(end)]);%ylim([0 1.5]);
    %text(ceil(x/2),0.4,num2str(length(clust)),'fontsize',25);hold on;
    title(['class ' num2str(i) '_' num2str(length(clust))],'fontsize',15,'interpreter','none')
    ylabel('△F/F0');xlabel('frame');
    ylim(y);xlim([frame.cs_start frame.per_cycle]);

    name=fieldnames(para_mean);
    txt_dis=string;
    txt_dis=[txt_dis;'number',': ',num2str(length(clust))];
    for ii=1:length(name)
        a = getfield(para_mean(i),name{ii});
        txt_dis=[txt_dis;name{ii},': ',num2str(a)];
    end
    txt = uicontrol('Style','text',...
    'Position',[10+5+(i-1)*100,50,100,220],...
    'fontsize',12,...
    'BackgroundColor','w',...
    'horizontalalignment','left',...
     'string',txt_dis);
end
hold off
y=[0 0.02];
figure,    
for i = 1:length(unique(cidx3))
    name=fieldnames(para_mean);
    txt_dis=string;
    txt_dis=[txt_dis;'number',': ',num2str(length(clust))];
    for ii=1:length(name)
        a = getfield(para_mean(i),name{ii});
        txt_dis=[txt_dis;name{ii},': ',num2str(a)];
    end
    txt = uicontrol('Style','text',...
    'Position',[10+5+(i-1)*100,50,100,200],...
    'fontsize',10,...
    'BackgroundColor','w',...
    'horizontalalignment','left',...
     'string',txt_dis);
end
hold off

for i = 1:length(unique(cidx3))
 fig{i}=show_spv_GUI(showspv_mean_class{i});
 seqwrite(showspv_mean_class{i},[actpath 'Class_' num2str(i)])
end



