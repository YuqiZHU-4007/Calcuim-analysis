%%%%%%%%%correlation
field=fieldnames(X4);
acti=[];indd1=[];B=[];indd2=[];
for ii=1:length(field)
    acti=[acti,getfield(X4,field{ii})];
    indd1=[indd1;[getfield(CS_acq_ind,field{ii})]];%Óã
    indd2=[indd2;ii*ones(size(getfield(X4,field{ii}),2),1)];%ÄÔÇø
    B=[B;getfield(XX1,field{ii})'];
end
B=B';
normmethod='zscore';

acti_norm=normalize(acti,normmethod);
figure,imagesc(acti_norm);

%ind=indd(:,1);

% [~,ind]=sort(indd1);
% acti=acti(:,ind);
% indd1=indd1(ind);
% indd2=indd2(ind);
% B=B(:,ind);
% 
% acti=X4.OB;
% B=XX1.OB;
% ind=CS_acq_ind.OB;
%acti_norm=normalize(acti,normmethod);


figure,Correlational_analysis(acti);
xticklabels({num2str(ind)}); %set(gca,'xtick',ind,'ytick',ind);

figure,Correlational_analysis(acti_norm);
h_gca.XTickLabel = [];h_gca.YTickLabel = [];

A=acti_norm';
%A=pca(acti);
%Classifcation
%ind_cut_trial=ind_cut_trial_total{2};
eucD = pdist(A,'Euclidean');
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

%%%%%%%%%%%%%%%%%%%k-means
classnum=3;
[cidx3,cmeans3,sumd3] = kmeans(A,classnum,'replicates',100,'display','final','Distance','sqEuclidean');%'cosine','sqEuclidean'
figure,[silh3,h] = silhouette(A,cidx3,'sqEuclidean');
figure,
x=1:size(B,1);
for i = 1:length(unique(cidx3))
    subplot(ceil(classnum/2),2,i)
    % subplot(1,5,i-5)
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
        'LineWidth',1,...
        'markersize',3);
    hold on
    plot(x,mean(B(:,clust),2),'-o','color',[1 0 0],...
        'LineWidth',2,...
        'markersize',3);
    hold on
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
    text(5,0.4,num2str(length(clust)),'fontsize',25);hold on;
    ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
    ylim([-0.05 0.2])
end
hold off

figure,
x=1:size(B,1);
for i = 1:length(unique(cidx3))
    %subplot(ceil(classnum/2),2,i)
    %subplot(2,ceil(classnum/2),i)
    subplot(4,1,i)
    clust = find(cidx3==i);
    length(clust)
    plot(A(clust,:)','b',...
        'LineWidth',1);
    hold on
    plot(mean(A(clust,:)),'r',...
        'LineWidth',3);
    hold on
    line([(8:8:size(A,2))',(8:8:size(A,2))'],[-0.05;0.05],'color','k','linewidth',2)
    set(gca,'FontSize',20);box off;
    ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
    ylim([-0.05 0.05])
end
hold off

T=[];T2=[];
figure,
color=[1 0 0; 0 1 0; 0.5 0.5 0.5; 0 0 1];
colorfiled=[0,0,0;...
    0.3,0.75,0.93;...
    0,0.45,0.74;...
    0,0,1;...
    0.93,0.69,0.13;...
    0.85,0.33,0.10;...
    1,0,0;...
    0.64,0.08,0.18;...
    0.49,0.18,0.56;...
    1,0,1;];
for i = 1:length(unique(cidx3))
    subplot(ceil(classnum/2),2,i)
    clust = find(cidx3==i);
    for jj=1:length(unique(indd2))
        ind_filed=find(indd2==jj);
        inter=intersect( clust, ind_filed);
        T(jj,i)=length(inter);
        if ~isempty(inter)
            plot(x,B(:,inter),'-o','color',colorfiled(jj,:,:));hold on;
        end
        for jjj=1:length(unique(indd1))
            ind_fish=find(indd1==jjj);
            inter_fish=intersect( inter, ind_fish);
            T2(jjj,i,jj)=length(inter_fish);
%             if ~isempty(inter_fish)
%                 plot(x,B(:,inter_fish),'-o','color',color(jjj,:,:));hold on;
%             end
        end
    end
end

T3={};
for ii=1:4
    T3{ii,1}(:,:)=T2(ii,:,:);
    T3{ii,1}=[T3{ii,1}]';
end


