%demo of network measures
%zyq,20190917
%refer to BCT
clc;
clear all;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath,'/Path.mat'])
load(fullfile('H:\1.Test US\2.Tail free¡ª¡ªData from 117\20220803\fish2','brain_region_related_statistic.mat'));
sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
group={'L','C','NL','FL'};
A=3;B=4;
thr=0.5;method_thr= 'absolute';'proportional';
load('H:\3.Juvenile reference brain\registration to templete\ÄÔÇø·Ö¸î\segmentation_file_0525_DSregion_mask.mat');
loc=[];
for ii=1:length(Label)
    a=getfield(reg_loc,strrep(Label{ii},' ','_'));
    loc(ii,:)=mean(a,1);
end
figure,
scatter3(loc(:,1),loc(:,2),loc(:,3))
%% load C
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
    C={}; C_W={};
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\',group{batchi},'\','avg']);
    path=Path{batchi};
    for fishi=1:length(path)
        for ii=1:10
            ind=find(squeeze(sum(sum(fc_region_all{batchi}(:,:,:,ii,fishi),1,'omitnan'),2,'omitnan'))~=0)
            a=squeeze(mean(fc_region_all{batchi}(:,:,ind,ii,fishi),3,'omitnan'));
            C{ii}{fishi}=a;
            switch method_thr
                case 'proportional'
                    a=threshold_proportional(abs(C{ii}{fishi}), thr);
                case  'absolute'
                    a=threshold_absolute(abs(C{ii}{fishi}), thr);%
            end
            a(isnan(a))=0; 
            C_W{ii}{fishi} = a;
        end
    end
    %% Degree and Similarity
    h=figure('Name','Degree Distribution','position',[1,1,1920,980]);
    Degree_region=[];
    for zz=1:length(C)
        counts=[];centers=[];
        for fishi=1:length(path)
            CIJ=C_W{zz}{fishi};%CIJ(isnan(CIJ))=0;
            [deg] = degrees_und(CIJ);
            a=CIJ;a(a~=0)=1;Degree_region(:,zz,fishi)=sum(a,1);
            [counts(:,fishi),centers(:,fishi),] = hist(deg,[1:30]);
        end
        subplot(A,B,zz),
        x=centers(:,1);y=counts;
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        bar(x,m);hold on;
        scatter(reshape(repmat(x,1,size(y,2)),1,[]),reshape(y,1,[]),10,'k','filled');hold on;
        errorbar(x,m,sem,'LineStyle','none','color','k' );hold on;
        xlabel('Node Degree');xlim([0 30]);
        title( sessionx{zz});ylim([0 10])
    end
    saveas(h,[p,'\','Degree Distribution'],'jpeg');
    %h=figure;plot(b);set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',45);ylabel('Degree');
    %saveas(h,[p,'\','Degree Distribution'],'jpeg');
    h=figure('Name','sort of degree','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),
        x=1:30;y=squeeze(Degree_region(:,zz,:));
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        [a,I]=sort(m,'descend');
        bar(x,m(I));hold on;
        scatter(reshape(repmat(x,1,size(y,2)),1,[]),reshape(y(I,:),1,[]),10,'k','filled');hold on;
        errorbar(x,m(I),sem(I),'LineStyle','none','color','k' );hold on;
        set(gca,'xtick',[1:30],'xticklabel',{Label{I}},'XTickLabelRotation',90)
        ylim([0 100]);
        xlabel('# Node');ylabel('Centrality of Edges');
        title( sessionx{zz});
    end
    saveas(h,[p,'\','sort of degree'],'jpeg');
    
    h=figure('position',[1,1,1920,980]);
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=Degree_region(ii,:,:);
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.5);hold on;
        plot(x,m,'r','linewidth',2);hold on;
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);title(Label{ii});ylim([0 30]);
    end
    saveas(h,[p,'\','Degree in region'],'jpeg');
    h=figure('position',[1,1,1920,980]);
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=(Degree_region(ii,:,:)-Degree_region(ii,1,:));
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.5);hold on;
        plot(x,m,'r','linewidth',2);hold on;
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);title(Label{ii});
        %ylim([0 30]);
    end
    saveas(h,[p,'\','Change in Degree in region'],'jpeg');
    %% Centrality
    Edge_Bet_Centra={};Bet_Centra=[];max_edge_centrality=[];max_node_centrality=[];coreness=[];
    for fishi=1:length(path)
        for zz=1:length(C)
            CIJ=C_W{zz}{fishi};
            WIJ = weight_conversion(CIJ, 'lengths');
            [Edge_Bet_Centra{zz}(:,:,fishi),~]=edge_betweenness_wei(WIJ);
            Bet_Centra(:,zz,fishi)=betweenness_wei(WIJ);
            [coreness(:,zz,fishi),~] = kcoreness_centrality_bu(CIJ);
        end
    end
    Bet_Centra_n=[];
    for fishi=1:length(path)
        a=squeeze(Bet_Centra(:,:,fishi));a=reshape(a,1,[]);
        Bet_Centra_n(:,:,fishi)=reshape(normalize(a,'range',[0 100]),size(Bet_Centra,1),[]);
    end
    h=figure('Name','edge_betweenness_centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),title( sessionx{zz});
        imagesc(mean(Edge_Bet_Centra{zz},3),[0 10]);colorbar;%[min(Edge_Bet_Centra{zz}(:)) max(Edge_Bet_Centra{zz}(:))/2]
        max_edge_centrality(zz,1)=min(find(Edge_Bet_Centra{zz}(:) == max(Edge_Bet_Centra{zz}(:))));
        set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
        set(gca,'ytick',[1:30],'yticklabel',Label,'yTickLabelRotation',0);title( sessionx{zz});
    end
    saveas(h,[p,'\','edge_betweenness_centrality'],'jpeg');
    
    h=figure('Name','nodal_betweenness centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),
        x=1:30;y=squeeze(Bet_Centra_n(:,zz,:));
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        bar(x,m);hold on;
        scatter(reshape(repmat(x,1,size(y,2)),1,[]),reshape(y,1,[]),10,'k','filled');hold on;
        errorbar(x,m,sem,'LineStyle','none','color','k' );hold on;
        xlabel('# Node');ylabel('Centrality of Nodes');
        title( sessionx{zz});ylim([0 100]);
        %max_node_centrality(zz,1:length(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz)))))=(find(Bet_Centra(:,zz)== max(Bet_Centra(:,zz))));
        set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
    end
    saveas(h,[p,'\','nodal_betweenness centrality'],'jpeg');
    h=figure('Name','sort of nodal_betweenness centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),
        x=1:30;y=squeeze(Bet_Centra_n(:,zz,:));
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        [a,I]=sort(m,'descend');
        bar(x,m(I));hold on;
        scatter(reshape(repmat(x,1,size(y,2)),1,[]),reshape(y(I,:),1,[]),10,'k','filled');hold on;
        errorbar(x,m(I),sem(I),'LineStyle','none','color','k' );hold on;
        set(gca,'xtick',[1:30],'xticklabel',{Label{I}},'XTickLabelRotation',90)
        ylim([0 100]);
        xlabel('# Node');ylabel('Centrality of Nodes');
        title( sessionx{zz});
    end
    saveas(h,[p,'\','sort of nodal_betweenness centrality'],'jpeg');
    
    h=figure('Name','sort of edge_betweenness centrality','position',[1,1,1920,980]);
    for zz=1:length(C)
        subplot(A,B,zz),
        x=1:30;y=squeeze(sum(Edge_Bet_Centra{zz}(:,:,:),1));
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        [a,I]=sort(m,'descend');
        bar(x,m(I));hold on;
        scatter(reshape(repmat(x,1,size(y,2)),1,[]),reshape(y(I,:),1,[]),10,'k','filled');hold on;
        errorbar(x,m(I),sem(I),'LineStyle','none','color','k' );hold on;
        set(gca,'xtick',[1:30],'xticklabel',{Label{I}},'XTickLabelRotation',90)
        ylim([0 100]);
        xlabel('# Node');ylabel('Centrality of Edges');
        title( sessionx{zz});
    end
    saveas(h,[p,'\','sort of edge_betweenness centrality'],'jpeg');

    
    h=figure('position',[1,1,1920,980]);
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=Bet_Centra_n(ii,:,:);
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.5);hold on;
        plot(x,m,'r','linewidth',2);hold on;title(Label{ii});ylim([0 60]);
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10])
    end
    saveas(h,[p,'\','centrality in region'],'jpeg');
    
    for zz=1:10
        aa=mean(Edge_Bet_Centra{zz},3);
        G = graph(aa ,Label,'upper','omitselfloops');
        h=figure('position',[1,1,1920,980]);
        c=hot(floor(max(mean(Bet_Centra_n(:,zz,:),3)+3)));%GetColormap('hot',100);
        colormap(c);pl=plot(G);colorbar;%,'EdgeLabel',G.Edges.Weight
        %p.NodeCData=Bet_Centra(:,zz);
        pl.NodeColor = c(fix(mean(Bet_Centra_n(:,zz,:),3)+1),:);
        pl.NodeFontSize=15;pl.MarkerSize=15;
        pl.XData=[loc(:,1)];pl.YData=[loc(:,2)];pl.ZData=[loc(:,3)];
        a=table2array(G.Edges(:,2));pl.LineWidth =7*a/max(a);title(sessionx{zz});
        saveas(h,[p,'\','fc display-',sessionx{zz}],'jpeg');
    end
    
    for zz=2:10
        aa=mean(Edge_Bet_Centra{zz}-Edge_Bet_Centra{1},3);
        G = graph(aa ,Label,'upper','omitselfloops');
        h=figure('position',[1,1,1920,980]);
        c=hot(floor(max(mean(Bet_Centra_n(:,zz,:),3)+3)));%GetColormap('hot',100);
        colormap(c);pl=plot(G);colorbar;%,'EdgeLabel',G.Edges.Weight
        %p.NodeCData=Bet_Centra(:,zz);
        pl.NodeColor = c(fix(mean(Bet_Centra_n(:,zz,:),3)+1),:);
        pl.NodeFontSize=15;pl.MarkerSize=15;
        pl.XData=[loc(:,1)];pl.YData=[loc(:,2)];pl.ZData=[loc(:,3)];
        a=table2array(G.Edges(:,2));
        c2=[0.5 0.5 0.5;0 0 1];pl.EdgeColor=repmat(c2(1,:),length(a),1);
        pl.EdgeColor(find(a > 0),:)=repmat(c2(2,:),length(find(a>0)),1);
        pl.LineWidth =abs(5*a/max(a))+min(abs(5*a/max(a)))+0.1;title(sessionx{zz});
        saveas(h,[p,'\','fc display-setdiff1',sessionx{zz}],'jpeg');
    end
    
    for zz=1:10
        y=squeeze(Bet_Centra_n(:,zz,:));
        m=squeeze(mean(y,2,'omitnan'));sd=squeeze(std(y,[],2,'omitnan'));sem=sd./sqrt(size(y,2));
        [~,I]=sort(m,'descend');
        regionid=I(1:10);
        aa=mean(Edge_Bet_Centra{zz}(regionid,regionid,:),3);G = graph(aa ,{Label{regionid}},'upper','omitselfloops');
        h=figure('position',[1,1,1920,980]);
        c=hot(floor(max(mean(Bet_Centra_n(regionid,zz,:),3)+3)));%GetColormap('hot',100);
        colormap(c);pl=plot(G);colorbar;%,'EdgeLabel',G.Edges.Weight
        %p.NodeCData=Bet_Centra(:,zz);
        pl.NodeColor = c(fix(mean(Bet_Centra_n(regionid,zz,:),3)+1),:);
        pl.NodeFontSize=15;pl.MarkerSize=15;
        a=table2array(G.Edges(:,2));pl.LineWidth =7*a/max(a);title(sessionx{zz});
        pl.XData=[loc(regionid,1)];pl.YData=[loc(regionid,2)];pl.ZData=[loc(regionid,3)];
        saveas(h,[p,'\','fc display-top10',sessionx{zz}],'jpeg');
    end
    save([p,'\fc.mat'],'Degree_region','Edge_Bet_Centra','Bet_Centra','Bet_Centra_n','coreness','C','C_W');
    
    N = size(C{1},1);gamma=1; % this parameter is an input for the computation of Eres.
    %% modularity
    C_W_m=[];
    for ss=1:10
        aa=[];
        for zz=1:size(C_W{ss},2)
            aa(:,:,zz)=C_W{ss}{zz};
        end
        C_W_m(:,:,ss)=mean(aa,3);
    end
    cut=[20,27,28];
    C_W_m(cut,:,:)=[];C_W_m(:,cut,:)=[];
    pl.NodeColor(find(Ci==1),:) =repmat(c(3,:),length(find(Ci==1)),1);
    pl.NodeColor(find(Ci==3),:) =repmat(c(1,:),length(find(Ci==3)),1);
    for zz=[1:10]
        aa=C_W{zz}{5};aa(cut,:)=[];aa(:,cut)=[];
        [Ci,Q]=modularity_und(aa,1);unique(Ci)
        aa=mean(Edge_Bet_Centra{zz},3);aa(cut,:)=[];aa(:,cut)=[];
        c=jet(length(unique(Ci)));
        G = graph(aa ,Label(setdiff(1:size(loc,1),cut)),'upper','omitselfloops');
        a=table2array(G.Edges(:,2));G = rmedge(G,find(a<0.1*max(a)));a=table2array(G.Edges(:,2));
        h=figure('position',[1,1,1920,980]);
        colormap(c);pl=plot(G);colorbar;
        pl.NodeColor =c(Ci,:);
        pl.LineWidth =7*a/max(a);title(sessionx{zz});
        pl.NodeFontSize=15;pl.MarkerSize=15;
        pl.XData=[loc(setdiff(1:size(loc,1),cut),1)];pl.YData=[loc(setdiff(1:size(loc,1),cut),2)];pl.ZData=[loc(setdiff(1:size(loc,1),cut),3)];
        %xlim([0 1024]);ylim([0 1024])
        set(gca,'visible','off')
        %saveas(h,[p,'\','modularity display-cut-',sessionx{zz}],'jpeg');
    end
%     regioni=[5 10 12 15]
%     for ii=1:10
%         aa=Edge_Bet_Centra{ii}(regioni,regioni);
% %         [X,Y,Z] = adjacency_plot_und(aa,loc(regioni,:));
% %         plot3(X,Y,Z)
%         G = graph(aa ,{Label{regioni}},'upper','omitselfloops');
%         h=figure('position',[1,1,1920,980]);
%         c=hot(floor(max(mean(Bet_Centra_n(regioni,zz,:),3)+3)));%GetColormap('hot',100);
%         colormap(c);pl=plot(G);colorbar;%,'EdgeLabel',G.Edges.Weight
%         %p.NodeCData=Bet_Centra(:,zz);
%         pl.NodeColor = c(fix(mean(Bet_Centra_n(regioni,zz,:),3)+1),:);
%         a=table2array(G.Edges(:,2));pl.LineWidth =7*a/max(a);title(sessionx{zz});
%         saveas(h,[p,'\','fc display-',sessionx{zz}],'jpeg');
%     end
%     
%     for zz=1:10
%         aa=mean(Edge_Bet_Centra{zz},3);
%         G = graph(aa ,Label,'upper','omitselfloops');
%         a=table2array(G.Edges(:,2));G = rmedge(G,find(a<0.3*max(a)));a=table2array(G.Edges(:,2));
%         h=figure('position',[1,1,1920,980]);        
%         c=hot(max(floor(mean(Bet_Centra_n(:,zz,:),3)))+1);%GetColormap('hot',100);
%         colormap(c);pl=plot(G);colorbar;
%         pl.NodeColor =c(floor(mean(Bet_Centra_n(:,zz,:),3)+1),:);
%         pl.LineWidth =7*a/max(a);title(sessionx{zz});
%         saveas(h,[p,'\','fc display-cut-',sessionx{zz}],'jpeg');
%     end
    close all;
end

clc;clear all;close all;     
Bet_Centra_n_all={};Degree_region_all={};
for batchi=1:4
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\',group{batchi},'\','avg']);
    load([p,'\fc.mat']);
    Bet_Centra_n_all{batchi}=Bet_Centra_n;
    Degree_region_all{batchi}=Degree_region;
end

path=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\']);
h=figure('position',[1,1,1920,980]);
c=[1 0 0;0.5 0.5 0.5;0 0 1;1 0 1];
for batchi=[4 2 1]
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=Bet_Centra_n_all{batchi}(ii,:,:);
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on;
        plot(x,m,'color',c(batchi,:),'linewidth',2);hold on;title(Label{ii});ylim([0 60]);
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);
       if batchi==1
        y1= squeeze(Bet_Centra_n_all{1}(ii,:,:));
        y2=squeeze( Bet_Centra_n_all{2}(ii,:,:));
        y3=squeeze( Bet_Centra_n_all{4}(ii,:,:));
        p=[];p2=[];
        for zz=1:10
            p(zz)=ranksum(y1(zz,:),y2(zz,:));
            p2(zz)=ranksum(y1(zz,:),y3(zz,:));
        end
        a=squeeze(mean(y1,3,'omitnan'));ind=find(p<0.05);
        scatter(ind,m(ind)+2*sem(ind),14,'r*');hold on;
        
        a=squeeze(mean(y3,3,'omitnan'));ind=find(p2<0.05);
        scatter(ind,zeros(length(ind),1),14,c(4,:),'*');hold on;
       end
    end
end
saveas(h,[path,'\','centrality in region_all'],'jpeg');

h=figure('position',[1,1,1920,980]);
for batchi=[4 2 1]
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\',group{batchi},'\','avg']);
    load([p,'\fc.mat']);
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=Degree_region_all{batchi}(ii,:,:);
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on;
        plot(x,m,'color',c(batchi,:),'linewidth',2);hold on;
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);title(Label{ii});ylim([0 30]);
        
        if batchi==1
            y1= squeeze(Degree_region_all{1}(ii,:,:));
            y2=squeeze( Degree_region_all{2}(ii,:,:));
            y3=squeeze( Degree_region_all{4}(ii,:,:));
            p=[];p2=[];
            for zz=1:10
                p(zz)=ranksum(y1(zz,:),y2(zz,:));
                 p2(zz)=ranksum(y1(zz,:),y3(zz,:));
            end
        a=squeeze(mean(y1,3,'omitnan'));ind=find(p<0.05);
        scatter(ind,m(ind)+2*sem(ind),14,'r*');hold on;
        
        a=squeeze(mean(y3,3,'omitnan'));ind=find(p2<0.05);
        scatter(ind,zeros(length(ind),1),14,c(4,:),'*');hold on;
        end
    end
end
saveas(h,[path,'\','Degree in region_all'],'jpeg');

h=figure('position',[1,1,1920,980]);
for batchi=[4 2 1]
    p=checkpath(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\fc\',method_thr,'\','thr-',num2str(thr),'\',group{batchi},'\','avg']);
    load([p,'\fc.mat']);
    for ii=1:30
        subplot(5,6,ii),
        x=1:10; y=(Degree_region_all{batchi}(ii,:,:)-Degree_region_all{batchi}(ii,1,:));
        m=squeeze(mean(y,3,'omitnan'));sd=squeeze(std(y,[],3,'omitnan'));sem=sd./sqrt(size(y,3));
        shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on;
        plot(x,m,'color',c(batchi,:),'linewidth',2);hold on;
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);title(Label{ii});
        
        if batchi==1
            y1= squeeze(Degree_region_all{1}(ii,:,:)-Degree_region_all{1}(ii,1,:));
            y2=squeeze(Degree_region_all{2}(ii,:,:)-Degree_region_all{2}(ii,1,:));
             y3=squeeze(Degree_region_all{4}(ii,:,:)-Degree_region_all{4}(ii,1,:));
            p=[];p2=[];
            for zz=1:10
                p(zz)=ranksum(y1(zz,:),y2(zz,:));
                p2(zz)=ranksum(y1(zz,:),y3(zz,:));
            end
        a=squeeze(mean(y1,3,'omitnan'));ind=find(p<0.05);
        scatter(ind,m(ind)+2*sem(ind),14,'r*');hold on;
        
        a=squeeze(mean(y3,3,'omitnan'));ind=find(p2<0.05);
        scatter(ind,zeros(length(ind),1),14,c(4,:),'*');hold on;
        end
    end
end
saveas(h,[path,'\','changes in Degree in region_all'],'jpeg');



