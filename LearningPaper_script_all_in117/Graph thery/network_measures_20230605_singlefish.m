%demo of network measures
%zyq,20190917
%refer to BCT
clc;
clear all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath,'/Path.mat'])
 load(fullfile('H:\1.Test US\2.Tail free¡ª¡ªData from 117\20220803\fish2','brain_region_related_statistic.mat'));
sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
group={'L','FL','NL','C'};
%% load C
C={}; C_W={};
A=3;B=4;
thr=0.4;method_thr='absolute';'proportional';
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
    path=Path{batchi};
    for fishi=1:length(path)
     nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\','thr-',num2str(thr),'\',group{batchi},'\',nn]);
    for ii=1:10
        ind=find(squeeze(sum(sum(fc_region_all{batchi}(:,:,:,ii,fishi),1,'omitnan'),2,'omitnan'))~=0)
        a=squeeze(mean(fc_region_all{batchi}(:,:,ind,ii,fishi),3,'omitnan'));
        C{ii}=a;
        switch method_thr
            case 'proportional'
                a=threshold_proportional(C{ii}, thr);
            case  'absolute'
                a=threshold_absolute(C{ii}, thr);%
        end
        a(isnan(a))=0; C_W{ii} = a;
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
        WIJ = weight_conversion(CIJ, 'lengths');
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
end
