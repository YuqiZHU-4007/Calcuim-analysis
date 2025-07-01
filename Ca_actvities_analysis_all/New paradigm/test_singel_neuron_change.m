load('')
num=xlsread('E:\A_Data_lightsheet\Data_vmat\20190223\fish3\20190223_fish3_ROI\location-20190223fish3','Raphe');
%[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
%num=[ones(1,48);1:48]';num=num(8:32,:);
num(find(num(:,4)==1),:)=[];
num_left=num(find(num(:,3)==1),:);
num_right=num(find(num(:,3)~=1),:);
ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
%%%%%%%%%correlation
acti=activities_preCS_dfdf_aftcorrect{1,1}(:,num(:,2));
acti_acq=[];
for ii=trial.acq(2):trial.acq(3)
    acti_acq=[acti_acq;acti((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.us_start-1,:)];
end
 figure,Correlational_analysis(acti_acq);
 set(gca,'xtick',1:size(acti_acq,2),'ytick',1:size(acti_acq,2));
 
a=num;acq_mean=[];normmethod='zscore';
for ii=1:size(a,1)
    acq_mean(:,ii)=area.CS_acq_block{a(ii,1)}(:,a(ii,2));
end
 acq_mean_norm=normlize_zyq_20190320(1,acq_mean,normmethod);
  acti_acq=normlize_zyq_20190320(1,acti_acq,normmethod);

 A=acti_acq;
 B=acq_mean;
figure,
x = categorical([ "Acq.-Bloc.1" "Acq.-Bloc.2" "Acq.-Bloc.3" "Acq.-Bloc.4" "Acq.-Bloc.5" ...
    "Acq.-Bloc.6" "Acq.-Bloc.7" "Acq.-Bloc.8" "Acq.-Bloc.9" "Acq.-Bloc.10"]','Ordinal',true);
x=1:size(B,1);
plot(x,A','-o',...
    'LineWidth',2,...
    'markersize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
ylabel('Integrated ¡÷F/F0');xlabel('Blocks');

%Classifcation
%ind_cut_trial=ind_cut_trial_total{2};
eucD = pdist(A','cosine');
clustTree = linkage(eucD,'average');
cophenet(clustTree,eucD)
figure,
[h,nodes] = dendrogram(clustTree,200);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
% [sum(ismember(nodes,[11 12 9 10])) sum(ismember(nodes,[6 7 8])) ...
%                   sum(ismember(nodes,[1 2 4 3])) sum(nodes==5)]
              
hidx = cluster(clustTree,'criterion','distance','cutoff',0.15,'depth',3);
unique(hidx)
for i = 1:length(unique(hidx))
    clust = find(hidx==i);
    figure,
    length(clust)
    idx=num(clust,:);acti=[];
    plot(x,B(:,clust),'-o',...
    'LineWidth',2,...
    'markersize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
    hold on
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
ylabel('Integrated ¡÷F/F0');xlabel('Blocks');

end
hold off

for i = 1:length(unique(hidx))
    clust = find(hidx==i);
    figure,
    length(clust)
    idx=num(clust,:);acti=[];
    %[area_CS_acq_block] = plot_meantrace_acq(mean(acti,2),fs,trial,frame);
    plot(x,mean(B(:,clust),2),'-o',...
    'LineWidth',2,...
    'markersize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
    hold on
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
end
hold off

for i = 1:length(unique(hidx))
    clust = find(hidx==i);
    length(clust)
    idx=num(clust,:);acti=[];
    figure,
    for jj=1:size(idx,1)
        acti(:,jj)=activities_preCS_dfdf_aftcorrect{idx(jj,1)}(:,idx(jj,2)); 
        acq_mean_cluster(:,jj)=area.CS_acq_block{idx(jj,1)}(:,idx(jj,2));
%         aa=xlsread(['E:\A_Data_lightsheet\Data_vmat\20190223\fish3\20190223_fish3_ROI\location\xlsx\','z',num2str(idx(jj,1),'%2d')]);
%         loca=aa(1,2+4*(idx(jj,2)-1));
    end
    %acti=mean(acti,2);
    acq_mean_cluster=mean(acq_mean_cluster,2);
    [~] = plot_meantrace_acq(mean(acti,2),fs,trial,frame,true,ind_cut_trial);
end
hold off

%%%%%%%%%%%%%%%%%%%k-means
[cidx3,cmeans3,sumd3] = kmeans(A',6,'replicates',6,'display','final','Distance','cosine');%'cosine','sqEuclidean'
figure,[silh3,h] = silhouette(A',cidx3,'cosine');
%figure,
for i = 1:length(unique(cidx3))
    clust = find(cidx3==i);
        length(clust)
    idx=num(clust,:);acti=[];
    figure,
    plot(x,B(:,clust),'-o',...
    'LineWidth',2,...
    'markersize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
    hold on
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
ylabel('Integrated ¡÷F/F0');xlabel('Blocks');

end
hold off
% for i = 1:length(unique(cidx3))
%     clust = find(cidx3==i);
%     figure,
%     length(clust)
%     idx=num(clust,:);acti=[];
%     %[area_CS_acq_block] = plot_meantrace_acq(mean(acti,2),fs,trial,frame);
%     plot(x,mean(B(:,clust),2),'-o',...
%     'LineWidth',2,...
%     'markersize',5,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r');
%     hold on
%     set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
% ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
% end
% hold off
for i = 1:length(unique(cidx3))
    clust = find(cidx3==i);
    figure,
    length(clust)
    idx=num(clust,:);acti=[];
    for jj=1:size(idx,1)
        acti(:,jj)=activities_preCS_dfdf_aftcorrect{idx(jj,1)}(:,idx(jj,2)); 
        acq_mean_cluster(:,jj)=area.CS_acq_block{idx(jj,1)}(:,idx(jj,2));
%         aa=xlsread(['E:\A_Data_lightsheet\Data_vmat\20190223\fish3\20190223_fish3_ROI\location\xlsx\','z',num2str(idx(jj,1),'%2d')]);
%         loca=aa(1,2+4*(idx(jj,2)-1));
    end
    %acti=mean(acti,2);
    acq_mean_cluster=mean(acq_mean_cluster,2);
    [~] = plot_meantrace_acq(mean(acti,2),fs,trial,frame,true,ind_cut_trial);
end
hold off


