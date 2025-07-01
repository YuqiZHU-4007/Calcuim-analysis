function plot_2D_node_clusnum   
C=graph.node_mean_act_training
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
        saveas(h,[input.savepath 'C_clusnum' num2str(num_clus(ii),'%02d'), '.eps'],savetype);
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
        saveas(h,[input.savepath 'Colorbar_clusnum' num2str(num_clus(ii),'%02d'), '.eps'],savetype);
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
    end
