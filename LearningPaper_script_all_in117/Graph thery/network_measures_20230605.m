%demo of network measures
%zyq,20190917
%refer to BCT
clc;
clear all;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath,'/Path.mat'])
 load(fullfile('H:\1.Test US\2.Tail free¡ª¡ªData from 117\20220803\fish2','brain_region_related_statistic.mat'));
sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
group={'L','FL','NL','C'};
load('H:\3.Juvenile reference brain\registration to templete\ÄÔÇø·Ö¸î\segmentation_file_0525_DSregion_mask.mat');
loc=[];
for ii=1:length(Label)
    a=getfield(reg_loc,strrep(Label{ii},' ','_'));
    loc(ii,:)=mean(a(ii,:),1);
end
%% load C
C={}; C_W={};
A=3;B=4;
thr=0.5;method_thr='absolute';'proportional';
fc_region_all={};
for batchi=1:4
    path=Path{batchi};
    for fishi=1:length(path)
        p=path{fishi};
        load([checkpath(fullfile(p,'fc')),'\fc_region.mat']);
        fc_region_all{batchi}(:,:,:,:,fishi)=fc_region;
    end
end
for batchi=1:4
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\',group{batchi}]);
    for ii=1:10
        ind=find(squeeze(sum(sum(fc_region_all{batchi}(:,:,:,ii,1),1,'omitnan'),2,'omitnan'))~=0)
        a=squeeze(mean(abs(fc_region_all{batchi}(:,:,ind,ii,:)),3,'omitnan'));
        C{ii}=squeeze(mean(a,3,'omitnan'));
        switch method_thr
            case 'proportional'
                a=threshold_proportional(abs(C{ii}), thr);
            case  'absolute'
                a=threshold_absolute(abs(C{ii}), thr);%
        end
        a(isnan(a))=0; 
        C_W{ii} = a;
    end
    h=figure('position',[1,1,1920,980]);
    for ii=1:10
    subplot(A,B,ii),imagesc(abs(C_W{ii}),[0,1]);colorbar;
    set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
    set(gca,'ytick',[1:30],'yticklabel',Label,'yTickLabelRotation',0);
    set(gca,'fontsize',10);title(sessionx{ii});
    end
    %close(h);            
    saveas(h,[p,'\','correlation map'],'jpeg');
    %saveas(h,[checkpath(fullfile(p,'fc')),'\',sessionx{ii}],'jpeg');
    N = size(C{1},1);lambda = 0.5; % this parameter is an input for the computation of Eres.
    %% Degree and Similarity
    h=figure('Name','Degree Distribution','position',[1,1,1920,980]);
    b=[];
    for zz=1:length(C)
        CIJ=C_W{zz};%CIJ(isnan(CIJ))=0;
        [deg] = degrees_und(CIJ);
        a=CIJ;a(a~=0)=1;b(:,zz)=sum(a,1);
        subplot(A,B,zz),hist(deg,[0:30]);xlim([0 30])
        xlabel('Node Degree');
        title( sessionx{zz});ylim([0 8])
    end
    saveas(h,[p,'\','Degree Distribution'],'jpeg');
%     h=figure;plot(b);set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',45);ylabel('Degree');
%     saveas(h,[p,'\','Degree Distribution'],'jpeg');
    h=figure('position',[1,1,1920,980]);
    for ii=1:30
        subplot(5,6,ii),
        plot(b(ii,:),'linewidth',2);title(Label{ii});ylim([0 30]);
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10])
    end
   saveas(h,[p,'\','Degree in region'],'jpeg');
    %% Centrality
    Edge_Bet_Centra={};Bet_Centra=[];max_edge_centrality=[];max_node_centrality=[];coreness=[];
    for zz=1:length(C)
        CIJ=C_W{zz};
        WIJ =weight_conversion(CIJ, 'lengths');
        %BC=betweenness_wei(WIJ);
        [Edge_Bet_Centra{zz},~]=edge_betweenness_wei(WIJ);
        Bet_Centra(:,zz)=betweenness_wei(WIJ);
        [coreness(:,zz),~] = kcoreness_centrality_bu(CIJ); 
    end
    h=figure('Name','edge_betweenness_centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),title( sessionx{zz});
        imagesc(Edge_Bet_Centra{zz},[0 10]);colorbar;%[min(Edge_Bet_Centra{zz}(:)) max(Edge_Bet_Centra{zz}(:))/2]
        max_edge_centrality(zz,1)=min(find(Edge_Bet_Centra{zz}(:) == max(Edge_Bet_Centra{zz}(:))));
        set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
        set(gca,'ytick',[1:30],'yticklabel',Label,'yTickLabelRotation',0);title( sessionx{zz});
    end
    saveas(h,[p,'\','edge_betweenness_centrality'],'jpeg');
    h=figure('Name','nodal_betweenness centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),bar(Bet_Centra(:,zz));ylim([0 100]);
        xlabel('# Node');ylabel('Centrality of Nodes');
        title( sessionx{zz});
        max_node_centrality(zz,1:length(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz)))))=(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz))));
        set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
    end
    saveas(h,[p,'\','nodal_betweenness centrality'],'jpeg');
    h=figure('Name','sort of nodal_betweenness centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),
        [a,I]=sort(Bet_Centra(:,zz),'descend');
        bar(a);set(gca,'xtick',[1:30],'xticklabel',{Label{I}},'XTickLabelRotation',90)
        ylim([0 100]);
        xlabel('# Node');ylabel('Centrality of Nodes');
        title( sessionx{zz});
        max_node_centrality(zz,1:length(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz)))))=(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz))));    
    end
    saveas(h,[p,'\','sort of nodal_betweenness centrality'],'jpeg');
    Label(max_node_centrality(:,1))
      
    h=figure('position',[1,1,1920,980]);
    for ii=1:30
        subplot(5,6,ii),
        plot(Bet_Centra(ii,:),'linewidth',2);title(Label{ii});ylim([0 100]);
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10])
    end
    saveas(h,[p,'\','centrality in region'],'jpeg');

    for zz=1:10
        G = graph( Edge_Bet_Centra{zz},Label,'upper','omitselfloops');
        h=figure('position',[1,1,1920,980]);
        c=hot(max(Bet_Centra(:,zz)+3));%GetColormap('hot',100);
        colormap(c);pl=plot(G);colorbar;%,'EdgeLabel',G.Edges.Weight
        %p.NodeCData=Bet_Centra(:,zz);
        pl.NodeColor = c(Bet_Centra(:,zz)+1,:);
        a=table2array(G.Edges(:,2));pl.LineWidth =7*a/max(a);title(sessionx{zz});
        saveas(h,[p,'\','fc display-',sessionx{zz}],'jpeg');
    end
    
    for zz=1:10
        G = graph( Edge_Bet_Centra{zz},Label,'upper','omitselfloops');
        a=table2array(G.Edges(:,2));G = rmedge(G,find(a<0.3*max(a)));a=table2array(G.Edges(:,2));
        h=figure('position',[1,1,1920,980]);        c=hot(max(Bet_Centra(:,zz)+3));%GetColormap('hot',100);
        colormap(c);pl=plot(G);colorbar;
        pl.NodeColor = c(Bet_Centra(:,zz)+1,:);        
        pl.LineWidth =7*a/max(a);title(sessionx{zz});
        saveas(h,[p,'\','fc display-cut-',sessionx{zz}],'jpeg');
    end
close all;
end



% h=figure('name','Hub and Edges');
% for zz=1:length(thr)
%     subplot(A,B,zz);
%     jj=max_node_centrality(zz,1);    ii=find(C{zz}(jj,:)~= 0);
%     scatter(input.nodes(:,2),input.nodes(:,1),20,[0.5 0.5 0.5],'filled',...
%         'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
%     line([input.nodes(repmat(jj,1,length(ii)),2),input.nodes(ii,2)]',[input.nodes(repmat(jj,1,length(ii)),1),input.nodes(ii,1)]',...
%         'linewidth',1,'color','r');hold on;
%     scatter(input.nodes(ii,2),input.nodes(ii,1),20,'k','filled',...
%         'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5,'linewidth',2);hold on;
%     scatter(input.nodes(jj,2),input.nodes(jj,1),40,'k','*',...
%         'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5,'linewidth',2);hold on;
%     [lix,liy]=ind2sub(size(C{zz}),max_edge_centrality(zz,1));
%     line([input.nodes(lix,2),input.nodes(liy,2)],[input.nodes(lix,1),input.nodes(liy,1)],...
%         'linewidth',2,'color','y');hold on;
%     axis equal;
%     ax=gca;
%     set(ax,'visible','off');
%     title(num2str(thr(zz)));
% end
%
% %% Paths and Distances
% SPL={};hops={};Pmat={};
% Dist_M={};Num_edge={};
% lambda=[];
% efficiency=[];
% ecc=[];
% radius=[];
% diameter=[];
% avg_SPL=[];
% for zz=1:length(thr)
%     CIJ=C{zz};
%     WIJ = weight_conversion(CIJ, 'lengths');
%     %[Dist_M{zz},Num_edge{zz}]=distance_wei(WIJ);
%     [SPL{zz},hops{zz},Pmat{zz}] = distance_wei_floyd(WIJ);%
%     [lambda(zz,1),efficiency(zz,1),ecc(:,zz),radius(zz,1),diameter(zz,1)] = charpath(SPL{zz});
% end
% h=figure('Name','Shortest Path length');
% for zz=1:length(thr)
%     y=SPL{zz};y(find(isinf(y)))=nan;
%     avg_SPL(zz,1)=mean(y(:),'omitnan');
%     subplot(A,B,zz),imagesc(y);colorbar
%     title(num2str(thr(zz)));
% end
% h=figure('Name','Number of edges');
% for zz=1:length(thr)
%     subplot(A,B,zz),imagesc(hops{zz});colorbar
%     title(num2str(thr(zz)));
% end
% %% Clustering and Community Structure
% clust_coef=[];transitivity=[];
% for zz=1:length(thr)
%     CIJ=C{zz};
%     WIJ = weight_conversion(CIJ, 'normalize');
%     clust_coef(:,zz)=clustering_coef_wu(WIJ); %W (all weights must be between 0 and 1)
%     transitivity(zz,:)=transitivity_wu(WIJ);
% end
% avg_clust_coef=mean(clust_coef)';
% % a=[min(W_nrm(:)) max(W_nrm(:))];
% % figure;imagesc(W_nrm(idx_ind,idx_ind),a);colormap('jet');colorbar;
% % avg_clustering_coef=mean(C);
% % Erout = rout_efficiency(W_nrm);
% % [Ci,Q]=modularity_und(W_nrm);
% % M=link_communities(W_nrm);
% %% Small world
% Small_wordness=[];
% avg_C=[];L=[];Small_wordness_MC=[];
% for zz=1:length(thr)
%     CIJ=C{zz};
%     CIJ = weight_conversion(CIJ, 'normalize');
%     deg = degrees_und(CIJ);  % degree distribution of undirected network
%     m = sum(deg)/2;
%     K = mean(deg); % mean degree of network
%     [expectedC,expectedL] = ER_Expected_L_C(K,size(CIJ,1));  % L_rand and C_rand
%     [Small_wordness(zz,1),avg_C(zz,1),L(zz,1)] = small_world_ness(CIJ,expectedL,expectedC,1);
%     [Lrand,CrandWS] = NullModel_L_C(size(CIJ,1),m,100,1);
%     Lrand_mean = mean(Lrand(Lrand < inf));
%     [Small_wordness_MC(zz,1),~,~] = small_world_ness(CIJ,Lrand_mean,mean(CrandWS),1);  % Using WS clustering coefficient
% end
% %% display
% %disp([prob_SPL,avg_clustering_coef,transitivity,Erout,Q])
% %table(thr',lambda,efficiency,radius,diameter,avg_clust_coef,avg_SPL,transitivity,Small_wordness,avg_C,L,Small_wordness_MC)
% table(thr',efficiency,transitivity,avg_C,L,Small_wordness_MC)
% save([gpath '\measures'],'thr','SPL','hops','Pmat','Dist_M','Num_edge','lambda','efficiency','ecc','radius','diameter',...
%     'avg_SPL','clust_coef','transitivity','avg_clust_coef','Edge_Bet_Centra','Bet_Centra','max_edge_centrality','max_node_centrality',...
%     'Small_wordness','avg_C','L','Small_wordness_MC');
% %% Density and Rentian Scaling
% [kden,N,K] = density_und(CIJ);
% %% Assortativity and Core Structure
% r = assortativity_wei(CIJ,flag);
% %% Efficiency and Diffusion
%
% %% Motifs and self-similarity