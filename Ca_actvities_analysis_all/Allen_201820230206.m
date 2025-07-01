clc;clear all;close all;
%% project activity into state-cue-behave axis
%refer to Allen_2018_science
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
for typei=1:4
    for batchi=1:length(Path{typei})
        p=char(Path{typei}(batchi));%'H:\1.Test US\2.Tail free——Data from 117\20210822\fish1\';
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'activities_aft_process.mat'));
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'))
        A=activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:);
        A_r=reshape(A,frame.per_cycle,trial.total,[]);
        %% 1:用population act构建基向量：平均每个condition下所有神经元活动得到N X 1矩阵；QR分解正交基 N X 3矩阵
        %state：post learning-pre learning CS前basline活动【1】；post learning-pre learning 整个trial活动【2】
        %Cue：post leanning CS期间活动 - pre learning CS期间活动【1】；post leanning nogo trial CS期间活动 - pre learning nogo trial CS期间活动
        %behavior： post leanning go trial - post leanning nogo trial CS期间活动
        go_trial={};nogo_trial={};trial_ind_session={};
        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
        for ii=1:trial.acq_block_num+2
            switch ii
                case 1
                    trial_ind_session{ii}=trial.hab(2):trial.hab(3);
                case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
                    trial_ind_session{ii}=(trial.hab(3)+(ii-2)*trial.acq_block_trial)+1 :(trial.hab(3)+(ii-1)*trial.acq_block_trial);
                case trial.acq_block_num+2
                    trial_ind_session{ii}=trial.test(2):trial.test(3);
            end
            
            go_trial{ii}=unique(CR(find(CR(:,1)<=trial_ind_session{ii}(end) & CR(:,1)>=trial_ind_session{ii}(1)),1));
            nogo_trial{ii}=setdiff(trial_ind_session{ii},go_trial{ii});
        end
        trial_ind=struct;period_ind=struct;Difference=struct;
        trial_ind.state{1}{1}=trial.hab(2):trial.hab(3);trial_ind.state{1}{2}=trial.test(2):trial.test(3);
        period_ind.state{1}=1:frame.per_cycle;
        period_ind.state{2}=1:frame.cs_start-1;
        trial_ind.cue{1}{1}=trial.hab(2):trial.hab(3);trial_ind.cue{1}{2}=trial.test(2):trial.test(3);
        trial_ind.cue{2}{1}=intersect(trial.hab(2):trial.hab(3),nogo_trial{1});trial_ind.cue{2}{2}=intersect(trial.test(2):trial.test(3),nogo_trial{end});
        period_ind.cue{1}=frame.cs_start:frame.us_start-1;
        trial_ind.behav{1}{1}=intersect(trial.test(2):trial.test(3),nogo_trial{end});trial_ind.behav{1}{2}=intersect(trial.test(2):trial.test(3),go_trial{end});
        period_ind.behav{1}=frame.cs_start:frame.us_start-1;
        % first subtracted the per-cell mean activity across all conditions (thirsty/sated/stim/washout).
        Difference.state{1}=abs(squeeze(mean(mean(A_r(period_ind.state{1},trial_ind.state{1}{2},:),2,'omitnan')-mean(A_r(period_ind.state{1},trial_ind.state{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.state{2}=abs(squeeze(mean(mean(A_r(period_ind.state{2},trial_ind.state{1}{2},:),2,'omitnan')-mean(A_r(period_ind.state{2},trial_ind.state{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.cue{1}=abs(squeeze(mean(mean(A_r(period_ind.cue{1},trial_ind.cue{1}{2},:),2,'omitnan')-mean(A_r(period_ind.cue{1},trial_ind.cue{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.cue{2}=abs(squeeze(mean(mean(A_r(period_ind.cue{1},trial_ind.cue{2}{2},:),2,'omitnan')-mean(A_r(period_ind.cue{1},trial_ind.cue{2}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.behav{1}=abs(squeeze(mean(mean(A_r(period_ind.behav{1},trial_ind.behav{1}{2},:),2,'omitnan')-mean(A_r(period_ind.behav{1},trial_ind.behav{1}{1},:),2,'omitnan'),1,'omitnan')));
        %% qr分解基向量
        W={};
        for ii=1:length(Difference.state)
            for jj=1:length(Difference.cue)
                for zz=1:length(Difference.behav)
                    avg_diff=cat(2,Difference.cue{jj},Difference.behav{zz},Difference.state{ii});
                    [W{ii,jj,zz},~] =qr(avg_diff,0);
                end
            end
        end
        %% 2: trial avg act across conditions: Pre/training/post x go/no-go
        X_avg_act=[];bin=2/fs.ca;
        for ii=1:length(trial_ind_session)
            for jj=1:2
                switch jj
                    case 1
                        a=go_trial{ii};
                    case 2
                        a=nogo_trial{ii};
                end
                b=trial_ind_session{ii};c=intersect(a,b);
                d=squeeze(mean(A_r(:,c,:),2,'omitnan'))';
                dd=reshape(d,size(d,1),bin,[]);
                X_avg_act(:,:,ii,jj)=squeeze(mean(dd,2,'omitnan'));
            end
        end
        %% 3：W’x计算projection到不同axis分量
        Projection_vect={};
        clr=hsv(length(trial_ind_session));
        figure('color','w'),kk=1;
        for ii=1;length(Difference.state)
            for jj=1;length(Difference.cue)
                for zz=1:length(Difference.behav)
                    %subplot(length(Difference.state),length(Difference.cue),kk),
                    s=[1,2,4,6,8,length(trial_ind_session)];line=[];
                    for cond_session=s
                        for cond_gonogo=1:2
                            switch cond_gonogo
                                case 1 %go
                                    type='-';
                                case 2  %no-go
                                   type='-'; %type='--';
                            end
                            Projection_vect{ii,jj,zz}(:,:,cond_session,cond_gonogo)= W{ii,jj,zz}'*squeeze(X_avg_act(:,:,cond_session,cond_gonogo));
                            x=Projection_vect{ii,jj,zz}(1,:,cond_session,cond_gonogo);
                            y= Projection_vect{ii,jj,zz}(2,:,cond_session,cond_gonogo);
                            z=Projection_vect{ii,jj,zz}(3,:,cond_session,cond_gonogo);
                            pp=plot3( x,y,z,'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;line(cond_session)=pp;
                            plot3( x(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),y(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),z(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;
                            scatter3(x(1),y(1),z(1),14,clr(cond_session,:),'filled');hold on;
                            n=3;scatter3(x(1:n:end),y(1:n:end),z(1:n:end),14,clr(cond_session,:),">",'filled');hold on;
                            U=x(2:end)-x(1:end-1);V=y(2:end)-y(1:end-1);w=z(2:end)-z(1:end-1);
                            b=0.01;
                            %                    quiver3(x(1:n:end-1),y(1:n:end-1),z(1:n:end-1),U(1:n:end)*b,V(1:n:end)*b,w(1:n:end)*b,0.1,'filled','color',clr(cond_session,:),'Markersize',6,'LineStyle',type,'MaxHeadSize',0.1,'Linewidth',2,'AutoScale','off','MarkerFaceColor',clr(cond_session,:));
                            hold on;
                        end
                    end
                    legend([line(s)],{num2str(s')},'Location','best')
                    xlabel('Cue');ylabel('Behav');zlabel('State');grid on;
                    %set(gca,'fontsize',16,'FontWeight','Bold')
                    kk=kk+1;
                end
            end
        end
        ylim([-1.5 1.5]);xlim([-2 2]);zlim([-1 1]);view([28,8]);grid off;title(nn);
        %% plot avg 分量
        save([p,'projection to StateCueMotion.mat'],'trial_ind','period_ind','Difference','W','X_avg_act','Projection_vect','-v7.3')
    end
end
for i = 1:72
   camorbit(-5,0,'camera')
drawnow
end
axis vis3d  %3维坐标系
for i=1:72
    camorbit(-5,0,'camera')
    drawnow
    M=getframe(gcf);
    nn=frame2im(M);
    [nn,cm]=rgb2ind(nn,256);
    if i==1
        imwrite(nn,cm,'out.gif','gif','LoopCount',inf,'DelayTime',0.1);%说明loopcount只是在i==1的时候才有用
    else
        imwrite(nn,cm,'out.gif','gif','WriteMode','append','DelayTime',0.1)%当i>=2的时候loopcount不起作用
    end
end
%% 不同脑区 不同type
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
Label={'L OB','R OB','L P','R P','L H','R H','L TeO','R TeO','L Np','R Np','L Otr','R Otr','TL','PO','Pr','PT','Th','rH','iH','cH','mPT','T','L TS','R TS','Va','CCe','L gT','R gT','R1','R2'};
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'}; 
projection_state_region={};projection_cue_region={};projection_behav_region={};
for typei=1:4
    for batchi=1:length(Path{typei})
        p=char(Path{typei}(batchi));%'H:\1.Test US\2.Tail free——Data from 117\20210822\fish1\';
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'para.mat'));load(fullfile(p,'projection to StateCueMotion.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'))
        kk=1;
        for regioni=1:length(Label)
            ind_in_regioni=find(strcmp(brain_region_id(:,2),Label(regioni)));
            s=[1:length(trial_ind_session)];
            % figure, %subplot(4,2,kk);
            for cond_session=s
                for cond_gonogo=1:2
                    switch cond_gonogo
                        case 1 %go
                            type='-';
                        case 2 %nogo
                            type='--';
                    end
                    avg_diff=cat(2,Difference.cue{1}(ind_in_regioni,:),Difference.behav{1}(ind_in_regioni,:),Difference.state{1}(ind_in_regioni,:));
                    w=qr(avg_diff,0);
                    x=w'*squeeze(X_avg_act(ind_in_regioni,:,cond_session,cond_gonogo));
                    projection_state_region{typei,batchi}(regioni,cond_session,cond_gonogo,:)=x(3,:);
                    projection_cue_region{typei,batchi}(regioni,cond_session,cond_gonogo,:)=x(1,:);
                    projection_behav_region{typei,batchi}(regioni,cond_session,cond_gonogo,:)=x(2,:);
                    %x=W{1,2,1}(ind_in_regioni,:)'*squeeze(X(ind_in_regioni,:,cond_session,cond_gonogo));
                    %             p=plot3( x(1,:),x(2,:),x(3,:),'linewidth',1.2,'color',clr(cond_session,:),'LineStyle',type); hold on;line(cond_session)=p;
                    %             plot3( x(1,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),x(2,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),x(3,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;
                    %             scatter3(x(1,1),x(2,1),x(3,1),14,clr(cond_session,:),'filled');hold on;
                    %             n=3;scatter3(x(1,1:n:end),x(2,1:n:end),x(3,1:n:end),14,clr(cond_session,:),">",'filled');hold on;
                end
            end
            %     legend([line(s)],[sessionx(s)],'Location','best')
            %     xlabel('Cue');ylabel('Behav');zlabel('State');grid on;title(Label(regioni));
            kk=kk+1;
        end
    end
end
projection_cue_all=[];projection_behav_all=[];projection_state_all=[];projection_cue_all_mean={};projection_behav_all_mean={};projection_state_all_mean={};
for typei=1
     s=[1:length(trial_ind_session)];
    for cond_session=s
        kk=1;
        for batchi=5;length(Path{typei})  
            NAME=projection_cue_region;
             b=squeeze( NAME{typei,batchi}(:,1,1,:));b(isnan(b))=0;
            bb=squeeze( NAME{typei,batchi}(:,1,2,:));bb(isnan(bb))=0;
            aa=squeeze( NAME{typei,batchi}(:,cond_session,1,:));aa(isnan(aa))=0;
            projection_cue_all(kk,cond_session,:,:)=squeeze(NAME{typei,batchi}(:,cond_session,2,:))- aa;
            NAME=projection_behav_region;
            b=squeeze( NAME{typei,batchi}(:,1,2,:));b(isnan(b))=0;
            bb=squeeze( NAME{typei,batchi}(:,1,1,:));bb(isnan(bb))=0;
            aa=squeeze( NAME{typei,batchi}(:,cond_session,2,:));aa(isnan(aa))=0;
            projection_behav_all(kk,cond_session,:,:)=(squeeze( NAME{typei,batchi}(:,cond_session,1,:))- aa);
            projection_state_all(kk,cond_session,:,:)=squeeze(projection_state_region{typei,batchi}(:,cond_session,2,:))- squeeze(projection_state_region{typei,batchi}(:,1,2,:))+squeeze(projection_state_region{typei,batchi}(:,cond_session,1,:))-squeeze( projection_state_region{typei,batchi}(:,1,1,:));
            %figure,imagesc(aa);
        kk=kk+1;
        end   
        projection_cue_all_mean{cond_session}=squeeze(mean( projection_cue_all(:,cond_session,:,:),1,'omitnan'));
        projection_behav_all_mean{cond_session}=squeeze(mean(projection_behav_all(:,cond_session,:,:),1,'omitnan'));
        projection_state_all_mean{cond_session}=squeeze(mean(projection_state_all(:,cond_session,:,:),1,'omitnan'));
    end
end
b=projection_state_all_mean;
a=[-0.3 0.3];
figure,subplot(1,4,1);imagesc(b{3},a);colormap(cool);colorbar;set(gca,'ytick',[1:30],'yticklabel',Label,'xticklabel',[],'fontsize',14)
subplot(1,4,2);imagesc(b{4},a);colormap(cool);colorbar;set(gca,'ytick',[1:30],'yticklabel',Label,'xticklabel',[],'fontsize',14)
subplot(1,4,3);imagesc(b{6},a);colormap(cool);colorbar;set(gca,'ytick',[1:30],'yticklabel',Label,'xticklabel',[],'fontsize',14)
%subplot(1,5,4);imagesc(projection_cue_all_mean{8},a);colormap(cool);colorbar;
subplot(1,4,4);imagesc(b{10},a);colormap(cool);colorbar;
set(gca,'ytick',[1:30],'yticklabel',Label,'xticklabel',[],'fontsize',14)
%% 
kk=1;
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'};
for regioni=1:length(Label)
    ind_in_regioni=find(strcmp(brain_region_id(:,2),Label(regioni)));
    s=[1,4,6,length(trial_ind_session)];x=[];
    for cond_session=s
        for cond_gonogo=1:2
            switch cond_gonogo
                case 1
                    type='--';
                case 2
                    type='-';
            end
            avg_diff=cat(2,Difference.cue{1}(ind_in_regioni,:),Difference.behav{1}(ind_in_regioni,:),Difference.state{1}(ind_in_regioni,:));
            w=qr(avg_diff,0);
            x(:,cond_session,1:3)=[w'*squeeze(X_avg_act(ind_in_regioni,:,cond_session,cond_gonogo))]';
        end 
    end
     figure,
     for kk=1:3
      subplot(1,3,kk);plot(squeeze(x(:,s,kk)));end
end