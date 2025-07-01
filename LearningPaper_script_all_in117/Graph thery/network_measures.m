%demo of network measures
%zyq,20190917
%refer to BCT
C=C_W;
lambda = 0.5; % this parameter is an input for the computation of Eres.
N = size(C{1},1);
EYE = logical(eye(N,N));
A=round(sqrt(length(thr)));B=round(length(thr)/sqrt(length(thr)));
%% weight_conversion
W_nrm={};
for zz=1:length(thr)
    CIJ=C{zz};
    W_nrm{zz} = weight_conversion(CIJ, 'normalize');
end
%% Degree and Similarity
figure('Name','Degree Distribution'),
for zz=1:length(thr)
    CIJ=C{zz};
    [deg] = degrees_und(CIJ);
    subplot(A,B,zz),hist(deg,[min(deg):max(deg)]);
    xlabel('Node Degree');
    xlim([1 max(deg)]);title(num2str(thr(zz)));
end
figure('Name','Strength');
for zz=1:length(thr)
    CIJ=C{zz};
    subplot(A,B,zz),
    [str] = strengths_und(CIJ);bar(str');
    xlabel('Strength')
    xlim([1 655]);title(num2str(thr(zz)));
end
%% Paths and Distances
SPL={};hops={};Pmat={};
Dist_M={};Num_edge={};
lambda=[];
efficiency=[];
ecc=[];
radius=[];
diameter=[];
avg_SPL=[];
for zz=1:length(thr)
    CIJ=C{zz};
    WIJ = weight_conversion(CIJ, 'lengths');
    %[Dist_M{zz},Num_edge{zz}]=distance_wei(WIJ);
    [SPL{zz},hops{zz},Pmat{zz}] = distance_wei_floyd(WIJ);%
    [lambda(zz,1),efficiency(zz,1),ecc(:,zz),radius(zz,1),diameter(zz,1)] = charpath(SPL{zz});
end
h=figure('Name','Shortest Path length');
for zz=1:length(thr)
    y=SPL{zz};y(find(isinf(y)))=nan;
    avg_SPL(zz,1)=mean(y(:),'omitnan');
    subplot(A,B,zz),imagesc(y);colorbar
    title(num2str(thr(zz)));
end
h=figure('Name','Number of edges');
for zz=1:length(thr)
    subplot(A,B,zz),imagesc(hops{zz});colorbar
    title(num2str(thr(zz)));
end
%% Clustering and Community Structure
clust_coef=[];transitivity=[];
for zz=1:length(thr)
    CIJ=C{zz};
    WIJ = weight_conversion(CIJ, 'normalize');
    clust_coef(:,zz)=clustering_coef_wu(WIJ); %W (all weights must be between 0 and 1)
    transitivity(zz,:)=transitivity_wu(WIJ);
end
avg_clust_coef=mean(clust_coef)';
% a=[min(W_nrm(:)) max(W_nrm(:))];
% figure;imagesc(W_nrm(idx_ind,idx_ind),a);colormap('jet');colorbar;
% avg_clustering_coef=mean(C);
% Erout = rout_efficiency(W_nrm);
% [Ci,Q]=modularity_und(W_nrm);
% M=link_communities(W_nrm);
%% Centrality
Edge_Bet_Centra={};Bet_Centra=[];max_edge_centrality=[];max_node_centrality=[];
for zz=1:length(thr)
    CIJ=C{zz};
    WIJ = weight_conversion(CIJ, 'lengths');
    %BC=betweenness_wei(WIJ);
    [Edge_Bet_Centra{zz},Bet_Centra(:,zz)]=edge_betweenness_wei(WIJ);
end
h=figure('Name','edge_betweenness_centrality');
for zz=1:length(thr)
    subplot(A,B,zz),
    imagesc(Edge_Bet_Centra{zz},[min(Edge_Bet_Centra{zz}(:)) max(Edge_Bet_Centra{zz}(:))/2]);colorbar
    title(num2str(thr(zz)));
    max_edge_centrality(zz,1)=min(find(Edge_Bet_Centra{zz}(:) == max(Edge_Bet_Centra{zz}(:))));
end
figure('Name','nodal_betweenness centrality'),
for zz=1:length(thr)
    subplot(A,B,zz),bar(Bet_Centra(:,zz));
    xlabel('# Node');title(num2str(thr(zz)));
    max_node_centrality(zz,1)=find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz)));
end
h=figure('name','Hub and Edges');
for zz=1:length(thr)
    subplot(A,B,zz);
    jj=max_node_centrality(zz,1);    ii=find(C{zz}(jj,:)~= 0);
    scatter(input.nodes(:,2),input.nodes(:,1),20,[0.5 0.5 0.5],'filled',...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
    line([input.nodes(repmat(jj,1,length(ii)),2),input.nodes(ii,2)]',[input.nodes(repmat(jj,1,length(ii)),1),input.nodes(ii,1)]',...
        'linewidth',1,'color','r');hold on;
    scatter(input.nodes(ii,2),input.nodes(ii,1),20,'k','filled',...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5,'linewidth',2);hold on;
    scatter(input.nodes(jj,2),input.nodes(jj,1),40,'k','*',...
        'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5,'linewidth',2);hold on;
    [lix,liy]=ind2sub(size(C{zz}),max_edge_centrality(zz,1));
    line([input.nodes(lix,2),input.nodes(liy,2)],[input.nodes(lix,1),input.nodes(liy,1)],...
        'linewidth',2,'color','y');hold on;
    axis equal;
    ax=gca;
    set(ax,'visible','off');    
    title(num2str(thr(zz)));
end
%% Small world
Small_wordness=[];
avg_C=[];L=[];Small_wordness_MC=[];
for zz=1:length(thr)
    CIJ=C{zz};
    CIJ = weight_conversion(CIJ, 'normalize');
    deg = degrees_und(CIJ);  % degree distribution of undirected network
    m = sum(deg)/2;
    K = mean(deg); % mean degree of network
    [expectedC,expectedL] = ER_Expected_L_C(K,size(CIJ,1));  % L_rand and C_rand
    [Small_wordness(zz,1),avg_C(zz,1),L(zz,1)] = small_world_ness(CIJ,expectedL,expectedC,1);
    [Lrand,CrandWS] = NullModel_L_C(size(CIJ,1),m,100,1);
    Lrand_mean = mean(Lrand(Lrand < inf));
    [Small_wordness_MC(zz,1),~,~] = small_world_ness(CIJ,Lrand_mean,mean(CrandWS),1);  % Using WS clustering coefficient    
end
%% display
%disp([prob_SPL,avg_clustering_coef,transitivity,Erout,Q])
%table(thr',lambda,efficiency,radius,diameter,avg_clust_coef,avg_SPL,transitivity,Small_wordness,avg_C,L,Small_wordness_MC)
table(thr',efficiency,transitivity,avg_C,L,Small_wordness_MC) 
save([gpath '\measures'],'thr','SPL','hops','Pmat','Dist_M','Num_edge','lambda','efficiency','ecc','radius','diameter',...
    'avg_SPL','clust_coef','transitivity','avg_clust_coef','Edge_Bet_Centra','Bet_Centra','max_edge_centrality','max_node_centrality',...
    'Small_wordness','avg_C','L','Small_wordness_MC');
%% Density and Rentian Scaling
[kden,N,K] = density_und(CIJ);
%% Assortativity and Core Structure
r = assortativity_wei(CIJ,flag);
%% Efficiency and Diffusion

%% Motifs and self-similarity

