%connection
%zyq-20190916
%undirect(Pearson Corr):refer to Richard,2018,BioRxiv
%direct(Granger Cau.):
addpath('F:\DUlab\FC analyse\Ca_actvities_analyze\Graph thery\BCT\2019_03_03_BCT');
[gname,gpath]=uigetfile('E:\A_Data_lightsheet\Data_huc\20190905\*.mat','graph');
load([gpath,gname]);
isspon2=true;
if isspon2 nbatch=2;else nbatch=1; end
savetype='fig';
isztrans=true;
thres_meth='absolute'; thr=0.5:0.05:0.9;thr_num=7;
input=graph;
outpath=checkpath([gpath 'graph_results']);
num_clus=2:10;
G=repmat([1:size(input.nodes,1)],size(input.nodes,1),1);
G2=repmat([1:size(input.nodes,1)]',1,size(input.nodes,1));
input = rmfield(input,'node_mean_act');
input = rmfield(input,'node_mean_act2');
input = rmfield(input,'node_mean_act_training');
%%
for batchi=1:nbatch
    batchi
    input.savepath=checkpath([outpath '\FC' num2str(batchi),'\']);
    switch batchi
        case 1
            input.node_mean_act=graph.node_mean_act;
        case 2
            input.node_mean_act=graph.node_mean_act2;
        case 3
            input.node_mean_act=graph.node_mean_act_training;
    end
    figure;
    n=1:size(input.node_mean_act,1);randi(size(input.node_mean_act,1),[1,1]);
    plot(input.node_mean_act(n,:)');xlim([1 size(input.node_mean_act,2)]);
    A=round(sqrt(length(thr)));B=round(length(thr)/sqrt(length(thr)));
    %% Pearson Corr
    C=corr(input.node_mean_act');
    Y =ones(size(C));Y=triu(Y);I=find(Y==0);I=[I;sub2ind([size(C)],1:size(C,1),1:size(C,1))'];
    a=[min(C(:)) max(C(:))];
    figure, imagesc(C,a);colormap('hot');colorbar;
    %% C V.S. distance
    X=squareform(pdist(input.nodes,'euclidean'));
    for ii=1
        %subplot(ceil(sqrt(length(thr))),ceil(sqrt(length(thr))),ii),
        Y=C;Y=triu(Y);
        X(I)=[];Y(I)=[];
        %     XX=sort(unique(X));YY=zeros(size(XX));kk=1;
        %    for jj=sort(unique(X));%min(X):max(X)
        %     YY(kk)=mean(Y(find(XX==jj)));
        %     kk=kk+1;
        %    end
        h=figure;scatter(X,Y,12,'b','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);hold on;
        %  plot(XX,YY,'color',[0.5 0.5 0.5],'linewidth',3);hold off;
        xlabel('Distance','fontsize',15);ylabel('Correlation','fontsize',15);
        line([0 max(X)],[0 0],'color','k');
        scatplot(X,Y,[],[],[],[],2,[]);hold on;
        saveas(h,[input.savepath 'C VS distance' '.',savetype],savetype);
    end
    %% subject by kmeans
    idx_ind_all=[];idx_all=[];
    for ii=9;
        [idx_all(:,ii),~,~,~] = kmeans(C',num_clus(ii),'Distance','correlation','Replicates',10);
        [~,idx_ind_all(:,ii)]=sort(idx_all(:,ii));
        %C=corr(input.node_mean_act(idx_ind,:)');
        h=figure('name','C');imagesc(C(idx_ind_all(:,ii),idx_ind_all(:,ii)),a);colormap('hot');colorbar;
        saveas(h,[input.savepath 'C_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        myColorMap = GetColormap('hsv_new',num_clus(ii));
        l=0;
        h=figure('name','Colorbar');
        for zz=1:num_clus(ii)
            ll=length(find(idx_all(:,ii)==zz));
            x = [l l+ll l+ll l];
            y = [0 0 1 1];
            l=l+ll;
            patch(x,y,myColorMap(zz,:));
        end
        xlim([1 length(idx_all(:,ii))]); ylim([0 1]);ax=gca;set(ax,'visible','off');
        saveas(h,[input.savepath 'Colorbar_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        h=figure('name','2D node');
        clr=myColorMap(idx_all(:,ii),:);
        temp=max(env.vol(:,:,1:24),[],3)';
        for zz=1;%:length(thr)
            %subplot(A,B,zz),
            imshow(temp,[min(temp(:)) max(temp(:))]);hold on;
            S=abs(mean(C,2))*100;
            scatter(input.nodes(:,2),input.nodes(:,1),S,clr,'filled',...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
            axis equal;
            ax=gca;
            set(ax,'visible','off');
        end
        saveas(h,[input.savepath '2D node_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        for zz=thr_num;
            %myColorMap = lines(size(input.node_mean_act(idx_ind,:),1));
            X=threshold_absolute(C, thr(zz));X(find(X>0))=1;
            h=figure('name',['circular_' num2str(thr(zz))]);
            circularGraph(X(idx_ind_all(:,ii),idx_ind_all(:,ii)),'Colormap',clr(idx_ind_all(:,ii),:));
            saveas(h,[input.savepath 'circular_clusnum' num2str(num_clus(ii),'%02d'),'_thr_', num2str(thr(zz)), '.',savetype],savetype);
        end
    end
    idx_ind=idx_ind_all(:,9);
    idx=idx_all(:,9);
    %% Granger
    %% Threshold
    C_T=C;
    switch thres_meth
        case 'proportional'
            if length(thr)>1
                h=figure('name','C_W');
                for ii=1:length(thr)
                    C_W{ii} = threshold_proportional(C_T, thr(ii));
                    subplot(A,B,ii);
                    a=[min(C_W{ii}(:)) max(C_W{ii}(:))];
                    imagesc(C_W{ii}(idx_ind,idx_ind),a);colormap('jet');colorbar;title(num2str(thr(ii)));
                end
                saveas(h,[input.savepath,'C_W', thres_meth ,'.',savetype],savetype);
            else
                C_W = threshold_proportional(C_T, thr);
                a=[min(C_W(:)) max(C_W(:))];
                h=figure('name','C_W');imagesc(C_W(idx_ind,idx_ind),a);colormap('jet');colorbar;title(num2str(thr));
                saveas(h,[input.savepath,'C_W', thres_meth ,'.',savetype],savetype);
            end
        case 'absolute'
            if length(thr)>1
                h=figure('name','C_W');
                for ii=1:length(thr)
                    C_W{ii} = threshold_absolute(C_T, thr(ii));
                    subplot(A,B,ii);
                    a=[min(C_W{ii}(:))  max(C_W{ii}(:))];
                    imagesc(C_W{ii}(idx_ind,idx_ind),a);colormap('jet');colorbar;title(num2str(thr(ii)));
                end
                saveas(h,[input.savepath,'C_W', thres_meth ,'.',savetype],savetype);
            else
                C_W = threshold_absolute(C_T, thr);
                a=[min(C_W(:)) max(C_W(:))];
                h=figure('name','C_W');imagesc(C_W(idx_ind,idx_ind),a);colormap('jet');colorbar;title(num2str(thr));
                saveas(h,[input.savepath,'C_W', thres_meth ,'.',savetype],savetype);
            end
    end
    %% Fisher transformation
    if isztrans
        C_T=C_W;
        if ~iscell(C_T)
            C_T=atan(C_T);%1/2.*log((1+C_T)./(1-C_T));
            a=[min(C_T(:)) max(C_T(:))];
            h=figure('name','C_T');imagesc(C_T(idx_ind,idx_ind),a);colormap('jet');colorbar;
                            saveas(h,[input.savepath,'C_T', thres_meth ,'.',savetype],savetype);
        else
            h=figure('name','C_T');
            for ii=1:length(thr)
                C_T{ii}=atan(C_T{ii});%1/2.*log((1+C_T{ii})./(1-C_T{ii}));
                a=[min(C_T{ii}(:))  max(C_T{ii}(:))];
                subplot(A,B,ii),
                imagesc(C_T{ii}(idx_ind,idx_ind),a);colormap('jet');colorbar;title(num2str(thr(ii)));
            end
                                                        saveas(h,[input.savepath,'C_T', thres_meth ,'.',savetype],savetype);
        end
    else
        C_T=C;
        a=[min(C_T(:)) max(C_T(:))];
        figure('name','C_T');imagesc(C_T(idx_ind,idx_ind),a);colormap('jet');colorbar;
    end
    %% Visualizaion
    idx_ind_all=[];idx_all=[];
    a=[0 max(C(:))];
    for ii=[4];
        [idx_all(:,ii),~,~,~] = kmeans(C',num_clus(ii),'Distance','correlation','Replicates',10);
%                 % change idx
%                 change_id=[1,2,5,4,3];
%                 for hh=1:length(unique(idx_all(:,ii)))
%                     aa{hh}=find(idx_all(:,ii)==hh);
%                 end
%                 for hh=1:length(unique(idx_all(:,ii)))
%                     idx_all(aa{hh},ii)=change_id(hh);
%                 end %change idx end
        [~,idx_ind_all(:,ii)]=sort(idx_all(:,ii));
        h=figure('name','C');imagesc(C(idx_ind_all(:,ii),idx_ind_all(:,ii)),a);colormap('hot');colorbar;
        saveas(h,[input.savepath 'C_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        myColorMap = [0.904,1,0.006;1,0.004,0.004;0.006,0.899,0;0.006,0.899,1;1,0.109,0.592];%GetColormap('hsv_new',num_clus(ii));[1 0 0;0 1 0];
        %[1,0.004,0.004;0.904,1,0.006;0.006,0.899,1;0.303,0.006,1;1,0.109,0.592]
        l=0;
        h=figure('name','Colorbar');
        for zz=1:num_clus(ii)
            ll=length(find(idx_all(:,ii)==zz));
            x = [l l+ll l+ll l];
            y = [0 0 1 1];
            l=l+ll;
            patch(x,y,myColorMap(zz,:));
        end
        xlim([1 length(idx_all(:,ii))]); ylim([0 1]);ax=gca;set(ax,'visible','off');
        saveas(h,[input.savepath 'Colorbar_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        h=figure('name','2D node');
        clr=myColorMap(idx_all(:,ii),:);
        temp=max(env.vol(:,:,1:24),[],3)';
        for zz=thr_num;%:length(thr)
            %subplot(A,B,zz),
            imshow(temp,[min(temp(:)) max(temp(:))]);hold on;
            S=abs(mean(C,2))*100;
            scatter(input.nodes(:,2),input.nodes(:,1),S,clr,'filled',...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
            %         ind=reshape(G,1,[]);ind2=reshape(G2,1,[]);
            %         lineS=abs(C_W{zz});lineS=reshape(lineS,1,[]);lineS(find(lineS==0))=nan;
            %         lineC=zeros(size(C_W{zz},1),size(C_W{zz},2),3);
            %         indC=find(C_W{zz}==0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0;lineC(indCx,indCy,2)=0;lineC(indCx,indCy,3)=0;
            %         indC=find(C_W{zz}>0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0.8;lineC(indCx,indCy,2)=0;lineC(indCx,indCy,3)=0.2;
            %         indC=find(C_W{zz}<0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0;lineC(indCx,indCy,2)=0.2;lineC(indCx,indCy,3)=0.8;
            %         lineC=reshape(lineC,[],1,3);
            %         for jj=1:length(ind)
            %             if isnan(lineS(jj))
            %                 continue;
            %             else
            %                 line([input.nodes(ind(jj),2) input.nodes(ind2(jj),2)],[input.nodes(ind(jj),1) input.nodes(ind2(jj),1)],'LineWidth',lineS(jj),'Color',lineC(jj,1,:));hold on;
            %             end
            %         end
            axis equal;
            ax=gca;
            set(ax,'visible','off');
        end
        saveas(h,[input.savepath '2D node_clusnum' num2str(num_clus(ii),'%02d'), '.',savetype],savetype);
        
        for zz=thr_num;
            h=figure('name',['circular_' num2str(thr(zz))]);
            %myColorMap = lines(size(input.node_mean_act(idx_ind,:),1));
            X=C_W{zz};X(find(X>0))=1;
            circularGraph(X(idx_ind_all(:,ii),idx_ind_all(:,ii)),'Colormap',clr(idx_ind_all(:,ii),:));
            saveas(h,[input.savepath 'circular_clusnum' num2str(num_clus(ii),'%02d'),'_thr_', num2str(thr(zz)), '.',savetype],savetype);
        end
    end
    %% 可不要    
     for jjj=1
%         %% 3-D connection
%         h=figure('name','3D connection');
%         for zz=1:length(thr)
%             subplot(A,B,zz),
%             X=C_T{zz};X(X>0)=1;
%             xlabel('X','fontsize',15);ylabel('Y','fontsize',15);zlabel('Z','fontsize',15);
%             for ii=1:size(C,1)-1
%                 for jj=ii+1:size(C,1)
%                     if X(ii,jj)==0
%                         break;
%                     else
%                         line([input.nodes(ii,1),input.nodes(jj,1)],[input.nodes(ii,2),input.nodes(jj,2)],[input.nodes(ii,3),input.nodes(jj,3)],...
%                             'linewidth',X(ii,jj),'color','r');hold on;
%                     end
%                 end
%             end
%             scatter3(input.nodes(:,1),input.nodes(:,2),input.nodes(:,3),15,'k','filled');hold on;
%             title(num2str(thr(zz)));
%         end
%         saveas(h,[input.savepath '\3D_connection.fig']);
%         h=figure('name','3D connection');
%         for zz=thr_num
%             X=C_T{zz};X(X>0)=1;
%             xlabel('X','fontsize',15);ylabel('Y','fontsize',15);zlabel('Z','fontsize',15);
%             for ii=1:size(C,1)-1
%                 for jj=ii+1:size(C,1)
%                     if X(ii,jj)==0
%                         break;
%                     else
%                         line([input.nodes(ii,1),input.nodes(jj,1)],[input.nodes(ii,2),input.nodes(jj,2)],[input.nodes(ii,3),input.nodes(jj,3)],...
%                             'linewidth',X(ii,jj),'color','r');hold on;
%                     end
%                 end
%             end
%             scatter3(input.nodes(:,1),input.nodes(:,2),input.nodes(:,3),15,'k','filled');hold on;
%             title(num2str(thr(zz)));
%             xlim([1 2048]);ylim([1 1024]);
%         end
%         saveas(h,[input.savepath '\3D_connection_thr_' num2str(thr(zz)) '.fig']);
%         %% 2-D plot
%         h=figure('name','2D connection');
%         myColorMap = GetColormap('hsv_new',length(unique(idx_ind)));
%         clr=myColorMap(idx,:);
%         temp=max(env.vol(:,:,1:24),[],3)';
%         for zz=6;%:length(thr)
%             %subplot(A,B,zz),
%             imshow(temp,[min(temp(:)) max(temp(:))]);hold on;
%             S=abs(mean(C,2))*100;
%             scatter(input.nodes(:,2),input.nodes(:,1),S,clr,'filled',...
%                 'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
%             %     ind=reshape(G,1,[]);ind2=reshape(G2,1,[]);
%             %     lineS=abs(C_W{zz});lineS=reshape(lineS,1,[]);lineS(find(lineS==0))=nan;
%             %     lineC=zeros(size(C_W{zz},1),size(C_W{zz},2),3);
%             %     indC=find(C_W{zz}==0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0;lineC(indCx,indCy,2)=0;lineC(indCx,indCy,3)=0;
%             %     indC=find(C_W{zz}>0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0.8;lineC(indCx,indCy,2)=0;lineC(indCx,indCy,3)=0.2;
%             %     indC=find(C_W{zz}<0);[indCx,indCy]=ind2sub(size(C_W{zz}),indC);lineC(indCx,indCy,1)=0;lineC(indCx,indCy,2)=0.2;lineC(indCx,indCy,3)=0.8;
%             %     lineC=reshape(lineC,[],1,3);
%             %     for ii=1:length(ind)
%             %         if isnan(lineS(ii))
%             %             continue;
%             %         else
%             %             line([input.nodes(ind(ii),2) input.nodes(ind2(ii),2)],[input.nodes(ind(ii),1) input.nodes(ind2(ii),1)],'LineWidth',lineS(ii),'Color',lineC(ii,1,:));hold on;
%             %         end
%             %     end
%             axis equal;
%             ax=gca;
%             set(ax,'visible','off');
%         end
%         %% circularGraph
%         % myLabel = cell(size(C,1));
%         % for i = 1:size(C,1)
%         %     myLabel{i} = num2str(i);
%         % end
%         for zz=thr_num;
%             %myColorMap = lines(size(input.node_mean_act(idx_ind,:),1));
%             X=C_W{zz};X(find(X>0))=1;
%             h=figure('name',['circular_' num2str(thr(zz))]);
%             circularGraph(X(idx_ind,idx_ind),'Colormap',clr(idx_ind,:));
%         end
    end
    %% 按脑区画
    h=figure('name','2D connection');
    myColorMap = GetColormap('hsv_new',length(unique(mask))+1);
    clr=myColorMap(input.node_label_region+1,:);
    temp=max(env.vol(:,:,1:24),[],3)';
    for zz=1;%:length(thr)
        %subplot(A,B,zz),
        imshow(temp,[min(temp(:)) max(temp(:))]);hold on;
        S=abs(mean(C,2))*100;
        scatter(input.nodes(:,2),input.nodes(:,1),S,clr,'filled',...
            'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.5);hold on;
        axis equal;
        ax=gca;
        set(ax,'visible','off');
    end
    [~,ii]=sort(input.node_label_region);   
    saveas(h,[input.savepath '2D connection_region', '.',savetype],savetype);
    for zz=thr_num;
        %myColorMap = lines(size(input.node_mean_act(idx_ind,:),1));
        X=C_W{zz};X(find(X>0))=1;
        h=figure('name',['circular_' num2str(thr(zz))]);
        circularGraph(X(ii,ii),'Colormap',clr(ii,:));
        saveas(h,[input.savepath '2D connection_region_thr_',num2str(thr(zz)), '.',savetype],savetype);
    end
    %% save
    save([input.savepath '\C.mat'],'C','C_T','C_W','idx_all','idx_ind_all');
end