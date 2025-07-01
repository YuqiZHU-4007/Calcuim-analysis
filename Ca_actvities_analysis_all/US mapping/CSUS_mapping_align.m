%% 时间尺度缩放因子 theta_trial
clc;clear all
global fs
global time
global frame
global frameb
global trial
global re_startpoint_sd
global startpoint
global y_3sd
savepath='H:\1.Test US\4.Hb activation\Hb_after_processing\';
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';

load([savepath '\Path']);
fix_trial=6;cont_win=2;alpha=0.05;align_trial=nan(2,length(Path{1}));
for batchi=[1,4]
    p=[];p2=[];p3=[];
    for fishi=1:length(Path{batchi})
        path=string(Path{batchi}(fishi));
        load(fullfile(path,'para.mat'));
        load(fullfile(path,'behavior_from_Results_of_alltheta.mat'));
        slidwin=6;trial_ind_session=[];
        for sessioni=1:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
            switch sessioni
                case 1
                    trial_ind_session(:,sessioni)=trial.hab(2):trial.hab(3);
                case mat2cell([2:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+1],1,ones(1,(trial.acq_block_num*trial.acq_block_trial-slidwin)+1));
                    trial_ind_session(:,sessioni)=trial.hab(3)+sessioni-1 :trial.hab(3)+sessioni+slidwin-2;
                case (trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
                    trial_ind_session(:,sessioni)=trial.test(2):trial.test(3);
            end
        end
        [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint_sd,[],[],[],1);
        x=cat(1,bef_cond.CS_shake,dur_cond.CS192_shake,aft_cond.CS_shake);
        y=p_value;
        p(size(trial_ind_session,2),fishi)=y.CS_15_6<=alpha & y.aftcond_befCS_CS_6<=alpha;
        p(2:size(trial_ind_session,2)-1,fishi)=y.dur_learnning_CS_6<=alpha & y.durcond_befCS_CS<=alpha;
        p2(size(trial_ind_session,2),fishi)= y.states_6<=alpha;
        p2(2:size(trial_ind_session,2)-1,fishi)= y.durcond_states_6<=alpha;
        for sessioni=1:size(trial_ind_session,2)
            p3(sessioni,fishi)=sum(x(trial_ind_session(:,sessioni)))./length(trial_ind_session(:,sessioni));
        end
        p_c=p;
        for cc=1:cont_win-1
            p_c(2+cont_win-1:size(trial_ind_session,2)-1,:)=p_c(2+cont_win-1:size(trial_ind_session,2)-1,:)+p([2+cont_win-1:size(trial_ind_session,2)-1]-cc,:);
        end
        ind=min(find(p3(2:size(trial_ind_session,2)-1,fishi)>=max(p3(2:size(trial_ind_session,2)-1,fishi)) & p(2:size(trial_ind_session,2)-1,fishi)>=1 ));% p_c(:,fishi)>=cont_win
        if ~isempty(ind) %
            if ind<fix_trial && p3(ind,fishi)>=0.5;
                align_trial(batchi,fishi)=ind;fix_trial;
            elseif ind>=fix_trial
                align_trial(batchi,fishi)=ind;
            end
        end
    end
    figure, imagesc(p3,[0 1]);text(1:length(Path{batchi}),align_trial(batchi,1:length(Path{batchi}))+1,'*');hold on;
end
save([savepath '\align_trial.mat'],'align_trial','slidwin','trial_ind_session');

%% 缩放activity
for batchi=[1 :4]
    for fishi=1:length(Path{batchi})
        if batchi==1 |batchi==4
            ind=align_trial(batchi,fishi)+slidwin-1;
        else
            ind=nan;
        end
        path=string(Path{batchi}(fishi))
        load(fullfile(path,'\activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(path,'\para.mat'));
        trial_hab=trial.hab(2):trial.hab(3);trial_test=trial.test(2):trial.test(3);
        go_trial_cs=unique(re_startpoint(find(re_startpoint(:,2)>=frameb.cs_start & re_startpoint(:,2)< frameb.us_start ),1))';
        nogo_trial_cs=setdiff(1:trial.total,go_trial_cs);
        go_trial_hab=intersect(go_trial_cs,[trial_hab])';go_trial_hab(size(go_trial_hab)+1:6,:)=nan;
        go_trial_tst=intersect(go_trial_cs,[trial_test])';go_trial_tst(size(go_trial_tst)+1:6,:)=nan;
        nogo_trial_hab=intersect(nogo_trial_cs,[trial_hab])';nogo_trial_hab(size(nogo_trial_hab)+1:6,:)=nan;
        nogo_trial_tst=intersect(nogo_trial_cs,[trial_test])'; nogo_trial_tst(size( nogo_trial_tst)+1:6,:)=nan;
        
        A_r=reshape(activities_preCS_dfdf_aftcorrect,frame.per_cycle,trial.total,[]);
        A_r= zscore(A_r,0,'all') ;
        if ind<fix_trial | ind>size(trial_ind_session,2)-1 | isnan(ind)
            align_win=reshape(trial.acq(2):trial.acq(3),[],fix_trial)
        elseif ind>=fix_trial & mod(ind,fix_trial)==0
            align_win=reshape(trial.acq(2):trial.acq(2)+ind-1,[],fix_trial)
        elseif ind>=fix_trial & mod(ind,fix_trial)~=0
            a=setdiff(trial.acq(2):trial.acq(2)+ind-1,ceil(linspace(trial.acq(2),trial.acq(2)+ind-2,mod(ind,fix_trial))));
            align_win=reshape(a,[],fix_trial)
        end
        align_win(size(align_win,1)+1:6,:)=nan;
        align_win_go=align_win;align_win_nogo=align_win;
        for ii=1:length(go_trial_cs)
            align_win_nogo(find(align_win==go_trial_cs(ii)))=nan;
        end
        align_win_go(find(~isnan(align_win_nogo)))=nan;
        activities_dfdf_align=nan(size(A_r,1),fix_trial,size(A_r,3));
        activities_dfdf_align_go=nan(size(A_r,1),fix_trial,size(A_r,3));
        activities_dfdf_align_nogo=nan(size(A_r,1),fix_trial,size(A_r,3));
        for ii=1:fix_trial
            a=align_win(:,ii);a(isnan(a))=[];
            activities_dfdf_align(:,ii,:)=squeeze(mean(A_r(:,a,:) ,2,'omitnan'));
            a=align_win_go(:,ii);a(isnan(a))=[];
            activities_dfdf_align_go(:,ii,:)=squeeze(mean(A_r(:,a,:) ,2,'omitnan'));
            a=align_win_nogo(:,ii);a(isnan(a))=[];
            activities_dfdf_align_nogo(:,ii,:)=squeeze(mean(A_r(:,a,:) ,2,'omitnan'));
        end
        
        a=mean(A_r(:,trial_hab,:) ,2,'omitnan');
        b=mean(A_r(:,trial_test,:) ,2,'omitnan');
        align_win=cat(2,trial_hab',align_win,trial_test');
        activities_dfdf_align=cat(2,a,activities_dfdf_align,b);
        a=mean(A_r(:,go_trial_hab(~isnan(go_trial_hab)),:) ,2,'omitnan');
        b=mean(A_r(:,go_trial_tst(~isnan(go_trial_tst)),:) ,2,'omitnan');
        activities_dfdf_align_go=cat(2,a,activities_dfdf_align_go,b);
        align_win_go=cat(2,go_trial_hab,align_win_go,go_trial_tst);
        a=mean(A_r(:,nogo_trial_hab(~isnan(nogo_trial_hab)),:) ,2,'omitnan');
        b=mean(A_r(:,nogo_trial_tst(~isnan(nogo_trial_tst)),:) ,2,'omitnan');
        activities_dfdf_align_nogo=cat(2,a,activities_dfdf_align_nogo,b);
        align_win_nogo=cat(2,nogo_trial_hab,align_win_nogo,nogo_trial_tst);
        save(fullfile(path, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo','activities_dfdf_align','activities_dfdf_align_go','activities_dfdf_align_nogo','-v7.3');
    end
end
%% cluster to superpixel
clc;clear all;close all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
for batchi=3:4
    for fishi=1:length(Path{batchi})
        path=string(Path{batchi}(fishi))
        load(fullfile(path, '/activities_dfdf_align.mat'));
        load(fullfile(path, '/env.mat'));
        A=cat(2,activities_dfdf_align,activities_dfdf_align_go,activities_dfdf_align_nogo);
        A(:,find(isnan(A(1,:,1))),:)=[];
        A = zscore(A,0,'all') ;
        A_r=reshape(A,size(A,1)*size(A,2),[]);
        %*******************************step1:cluster by spatial location**********
        res=[0.66,0.66,10];
        X=env.supervoxel(:,1:3);
        X(:,1)=X(:,1)*res(1);X(:,2)=X(:,2)*res(2);X(:,3)=(X(:,3)-1)*res(3);
        n=ceil(size(A,3)/10);
        clrmap = GetColormap('hsv_new',n);
        [idx_spacial,Center_spacial,sumd,D] = kmeans(X,n,'Distance','sqeuclidean','Replicates',10);
        within_clus=sqrt(sumd);m(1,1)=mean(within_clus);m(1,2)=std(within_clus);m(1,3)=length(within_clus);
        inter_clus=pdist(Center_spacial,'euclidean')';m(1,4)=mean(inter_clus);m(1,5)=std(inter_clus);m(1,6)=length(inter_clus);
%         clr=clrmap(idx_spacial,:);
%         temp = pdist(Center_spacial,'euclidean');Dist = squareform(temp);
%         figure;
%         scatter3(X(:,1),X(:,2),X(:,3),12,clr,'filled');hold on
%         plot3(Center_spacial(:,1),Center_spacial(:,2),Center_spacial(:,3),'kx','MarkerSize',15,'LineWidth',3);
%         xlabel('X');ylabel('Y');zlabel('Z');axis equal;
%         title 'Cluster Assignments and Centroids';hold off;
        h=figure;hist(idx_spacial,1:n);%xlim([1 n]);
        saveas(h,strcat(path,'\hist of cluster in step1.jpg'));    
        %******************************step2:cluster by activity*****************
        X=A_r;act_center=nan(length(unique(idx_spacial)),size(A_r,1));
        [coeff,score,latent,~,~,~] = pca(A_r); %A_r行观测值，列变量
        explained_var = cumsum(latent)/sum(latent);nPC=min(find(explained_var>0.9));
        n=nPC;
        for ii=unique(idx_spacial)'
            id=find(idx_spacial==ii);
            act_center(ii,:)=mean(X(:,id),2);
        end
        Cor=corr(act_center');
        [idx_functional,~,~,~] = kmeans(Cor',n,'Distance','correlation','Replicates',100);
        h=figure;hist(idx_functional,[1:n]);%xlim([1 n]);
        saveas(h,strcat(path, '\hist of cluster in step2.jpg'));
        h=figure;silhouette(Cor',idx_functional,'correlation');
        saveas(h,strcat(path, '\silhouette of cluster in step2.jpg'));
        [idxAs,idxAa_ind]=sort(idx_functional);
        h=figure;imagesc(Cor(idxAa_ind,idxAa_ind),[min(Cor(:)) max(Cor(:))]);colormap('hot');colorbar;
        saveas(h,strcat(path, '\corr of cluster in step2.jpg'));
        %***************************step3: merge and divide***********************
        idx_functional_adj=idx_functional;
        thr_divide=50;thr_merge=120;divided=ones(length(unique(idx_functional_adj')),1);a=[];
        for rep=1
            divided_idx=[];ii=0;
            while ii<=max(idx_functional_adj)-1;%unique(idxAadj')
                ii=ii+1;
                id=find(idx_functional_adj==ii);n=1;a(ii,1)=length(id);
                [~,Cn,sumd,~] = kmeans(Center_spacial(id,:),1,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
                if max(sumd)>=thr_divide
                    [idxn,Cn,sumd,~] = kmeans(Center_spacial(id,:),2,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
                    [counts,~]=hist(idxn,[1:max(idxn)]);
                    d_inter_clusters=squareform(pdist(Cn,'euclidean'));
                    [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
                    adj_merge=find(~isequal(merge_idx,merge_idy));
                    while  max(sumd)>=thr_divide && isempty(adj_merge) && min(counts)>=2; %mad>=thr1 && length(id)>=3 && squareform(pdist(C(id,:),'euclidean'))<=thr2
                        n=n+1;
                        if n>length(id)
                            warning([num2str(ii) 'not divided well!']);break;
                        else
                            [idxn,Cn,~,~] = kmeans(Center_spacial(id,:),n,'Distance','sqeuclidean','Replicates',10);
                            d_inter_clusters=squareform(pdist(Cn,'euclidean'));
                            [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
                            if ~isempty(merge_idx)
                                for jj=1:length(merge_idx)
                                    id1=find(idxn==merge_idx(jj));id2=find(idxn==merge_idy(jj));
                                    if ~isequal(id1,id2)
                                        idxn([id1,id2])=min(merge_idx(jj),merge_idy(jj));
                                    else
                                        break;
                                    end
                                end
                            end
                            [Cn,sumd] = FindCentroid_Direct(idxn,Center_spacial(id,:));sumd=sqrt(sumd);
                            d_inter_clusters=squareform(pdist(Cn,'euclidean'));
                            [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
                            adj_merge=find(~isequal(merge_idx,merge_idy));
                            [counts,~]=hist(idxn,[1:max(idxn)]);
                        end
                    end %while
                    divided_idx=id;%union(divided_idx,id);
                    adj=setdiff(1:length(idx_functional),union(id,divided_idx));
                    adj2=find(idx_functional_adj(adj)> ii);
                    %min(idxAadj(adj(adj2)))
                    idx_functional_adj(adj(adj2))=idx_functional_adj(adj(adj2))+max(idxn)-1;
                    %min(idxAadj(adj(adj2)))
                    idx_functional_adj(id)=idxn+idx_functional_adj(id)-1;
                    divided(ii,:)=n;
                    if length(unique(idx_functional_adj))~=length([1:max(idx_functional_adj)]')
                        disp true1
                        disp (num2str(ii))
                        break;
                    end
                else
                    divided_idx=id;%union(divided_idx,id);
                    divided(ii,:)=n;
                    %continue;
                    if length(unique(idx_functional_adj))~=length([1:max(idx_functional_adj)]')
                        disp false1
                        disp (num2str(ii));
                        break;
                    end
                end
            end
        end
        h=figure;
        subplot(2,1,1);hist(idx_functional,1:max(idx_functional));%xlim([1 max(idx_functional)]); max(idx_functional)
        subplot(2,1,2);hist(idx_functional_adj,1:max(idx_functional_adj));%xlim([1 max(idx_functional_adj)]);max(idx_functional_adj)
        saveas(h,strcat(path ,'\hist of cluster in step3.jpg'));
        [Center_final,~] = FindCentroid_Direct(idx_functional_adj,Center_spacial);
        h=figure;
        clrmap = GetColormap('hsv_new',length(unique(idx_functional_adj)));clr=clrmap(idx_functional_adj,:);
        X=Center_spacial;
        plot3(Center_final(:,1),Center_final(:,2),Center_final(:,3),'k*','MarkerSize',5,'LineWidth',1.5);
        xlabel('X');ylabel('Y');zlabel('Z');axis equal;
        title 'Cluster Assignments and Centroids';hold off;
        saveas(h,strcat(path, '\Centroids of cluster in step3.jpg'));
       %% check
        close all;
       %% save
       path=string(Path{batchi}(fishi))
        graph.nodes=Center_final;
        graph.node_mean_act=[];
        for ii=1:size(graph.nodes,1)
            graph.node_mean_act(ii,:)=mean(act_center(idx_functional_adj==ii,:),1);
        end
        graph.node_label=zeros(length(idx_spacial),1);
        for ii=1:length(idx_functional_adj)
            graph.node_label(find(idx_spacial==ii))=idx_functional_adj(ii);
        end
        graph.mean_act_catall=[];
        for ii=1:size(graph.nodes,1)
            graph.mean_act_catall(:,:,ii)=mean( A(:,:, graph.node_label==ii,:),3);
        end
        graph.idx_spacial=idx_spacial;
        graph.Center_spacial=Center_spacial;
        graph.idx_functional=idx_functional;
        graph.idx_functional_adj=idx_functional_adj;
        save(fullfile(path, '/activities_dfdf_align.mat'),'graph','-append','-v7.3');
    end
end

