clc;clear all;close all;
%% project activity into state-cue-behave axis
%refer to Allen_2018_science
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        p=char(Path{batchi}(fishi));%'H:\1.Test US\2.Tail free——Data from 117\20210822\fish1\';
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'))
        load(fullfile(p,'activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo')
        A=zscore(activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:),[],'all');
        A_r=reshape(A,frame.per_cycle,trial.total,[]);
        %% 1:用population act构建基向量：平均每个condition下所有神经元活动得到N X 1矩阵；QR分解正交基 N X 3矩阵
        %state：post learning-pre learning CS前basline活动【1】；post learning-pre learning 整个trial活动【2】
        %Cue：post leanning CS期间活动 - pre learning CS期间活动【1】；post leanning nogo trial CS期间活动 - pre learning nogo trial CS期间活动
        %behavior： post leanning go trial - post leanning nogo trial CS期间活动
        go_trial={};nogo_trial={};trial_ind_session={};
        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
        for ii=1:size(align_win,2)
            a=align_win(:,ii);a(isnan(a))=[];trial_ind_session{ii}=a;
            go_trial{ii}=unique(CR(find(CR(:,1)<=trial_ind_session{ii}(end) & CR(:,1)>=trial_ind_session{ii}(1)),1));
            nogo_trial{ii}=setdiff(trial_ind_session{ii},go_trial{ii});
        end
        trial_ind=struct;period_ind=struct;Difference=struct;
        %trial_ind.state{1}{1}=trial.hab(2):trial.hab(3);trial_ind.state{1}{2}=trial.test(2):trial.test(3);
        trial_ind.state{1}{2}=trial_ind_session{size(align_win,2)-1};trial_ind.state{1}{1}=trial_ind_session{2};
        period_ind.state{1}=1:frame.per_cycle;
        period_ind.state{2}=frame.us_start+4/fs.ca:frame.per_cycle;
        trial_ind.cue{1}{1}=trial.hab(2):trial.hab(3);trial_ind.cue{1}{2}=trial.test(2):trial.test(3);
        trial_ind.cue{2}{1}=intersect(trial.hab(2):trial.hab(3),nogo_trial{1});trial_ind.cue{2}{2}=intersect(trial.test(2):trial.test(3),nogo_trial{end});
        period_ind.cue{1}=frame.cs_start:frame.us_start-1;
        trial_ind.behav{1}{1}=intersect(trial.test(2):trial.test(3),nogo_trial{end});trial_ind.behav{1}{2}=intersect(trial.test(2):trial.test(3),go_trial{end});
        trial_ind.behav{2}{1}=intersect(trial.hab(2):trial.hab(3),nogo_trial{1});trial_ind.behav{2}{2}=intersect(trial.test(2):trial.test(3),go_trial{end});
        period_ind.behav{1}=frame.cs_start:frame.us_start-1;
        % first subtracted the per-cell mean activity across all conditions (thirsty/sated/stim/washout).
        Difference.state{1}=(squeeze(mean(mean(A_r(period_ind.state{1},trial_ind.state{1}{2},:),2,'omitnan')-mean(A_r(period_ind.state{1},trial_ind.state{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.state{2}=(squeeze(mean(mean(A_r(period_ind.state{2},trial_ind.state{1}{2},:),2,'omitnan')-mean(A_r(period_ind.state{2},trial_ind.state{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.cue{1}=(squeeze(mean(mean(A_r(period_ind.cue{1},trial_ind.cue{1}{2},:),2,'omitnan')-mean(A_r(period_ind.cue{1},trial_ind.cue{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.cue{2}=(squeeze(mean(mean(A_r(period_ind.cue{1},trial_ind.cue{2}{2},:),2,'omitnan')-mean(A_r(period_ind.cue{1},trial_ind.cue{2}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.behav{1}=(squeeze(mean(mean(A_r(period_ind.behav{1},trial_ind.behav{1}{2},:),2,'omitnan')-mean(A_r(period_ind.behav{1},trial_ind.behav{1}{1},:),2,'omitnan'),1,'omitnan')));
        Difference.behav{2}=(squeeze(mean(mean(A_r(period_ind.behav{1},trial_ind.behav{2}{2},:),2,'omitnan')-mean(A_r(period_ind.behav{1},trial_ind.behav{2}{1},:),2,'omitnan'),1,'omitnan')));
        
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
            for jj=1:3
                switch jj
                    case 1
                        a=go_trial{ii};
                    case 2
                        a=nogo_trial{ii};
                    case 3
                        a=trial_ind_session{ii};
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
        %figure('color','w'),kk=1;
        for ii=1:length(Difference.state)
            for jj=1:length(Difference.cue)
                for zz=1:length(Difference.behav)
                    %subplot(length(Difference.state),length(Difference.cue),kk),
                    s=[1:length(trial_ind_session)];line=[];
                    for cond_session=s
                        for cond_gonogo=1:3
                            switch cond_gonogo
                                case 1 %go
                                    type='-';
                                case 2  %no-go
                                    type='-'; %type='--';
                                case 3  %all
                                    type='-'; %type='--';
                            end
                            Projection_vect{ii,jj,zz}(:,:,cond_session,cond_gonogo)= W{ii,jj,zz}'*squeeze(X_avg_act(:,:,cond_session,cond_gonogo));
                            %                             x=Projection_vect{ii,jj,zz}(1,:,cond_session,cond_gonogo);x=normalize(x,"range");
                            %                             y= Projection_vect{ii,jj,zz}(2,:,cond_session,cond_gonogo);y=normalize(y,"range");
                            %                             z=Projection_vect{ii,jj,zz}(3,:,cond_session,cond_gonogo);z=normalize(z,"range");
                            %                             pp=plot3( x,y,z,'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;line(cond_session)=pp;
                            %                             plot3( x(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),y(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),z(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;
                            %                             scatter3(x(1),y(1),z(1),14,clr(cond_session,:),'filled');hold on;
                            %                             n=3;scatter3(x(1:n:end),y(1:n:end),z(1:n:end),14,clr(cond_session,:),">",'filled');hold on;
                            %                             U=x(2:end)-x(1:end-1);V=y(2:end)-y(1:end-1);w=z(2:end)-z(1:end-1);
                            %                             b=0.01;
                            %                             %                    quiver3(x(1:n:end-1),y(1:n:end-1),z(1:n:end-1),U(1:n:end)*b,V(1:n:end)*b,w(1:n:end)*b,0.1,'filled','color',clr(cond_session,:),'Markersize',6,'LineStyle',type,'MaxHeadSize',0.1,'Linewidth',2,'AutoScale','off','MarkerFaceColor',clr(cond_session,:));
                            %                             hold on;
                        end
                    end
                    %                     legend([line(s)],{num2str(s')},'Location','best')
                    %                     xlabel('Cue');ylabel('Behav');zlabel('State');grid on;
                    %                     %set(gca,'fontsize',16,'FontWeight','Bold')
                    %                     kk=kk+1;
                end
            end
        end
        %  ylim([-1 1]);xlim([-1 1]);zlim([-1 1]);view([28,8]);grid off;title(nn);
        %% plot avg 分量
        save([p,'projection to StateCueMotion_align.mat'],'bin','trial_ind','period_ind','Difference','W','X_avg_act','Projection_vect','trial_ind_session','go_trial','nogo_trial','-v7.3')
    end
end
projection_state_all={};projection_cue_all={};projection_behav_all={};
projection_state_all_mean={};projection_cue_all_mean={};projection_behav_all_mean={};
ii=2;jj=2;zz=2;
for cond_gonogo=1:3;
    for batchi=1:4
        for fishi=num{batchi}
            p=char(Path{batchi}(fishi));
            load([p,'projection to StateCueMotion_align.mat'],'Projection_vect');
            a=squeeze(Projection_vect{ii,jj,zz}(1,:,:,cond_gonogo));a=preprocessfunc(a,2);projection_cue_all{batchi,cond_gonogo}(:,:,fishi)= a;
            a=squeeze(Projection_vect{ii,jj,zz}(2,:,:,cond_gonogo));a=preprocessfunc(a,2);projection_behav_all{batchi,cond_gonogo}(:,:,fishi)= a;
            a=squeeze(Projection_vect{ii,jj,zz}(3,:,:,cond_gonogo));a=preprocessfunc(a,2);projection_state_all{batchi,cond_gonogo}(:,:,fishi)=a;
        end
        projection_cue_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_cue_all{batchi,cond_gonogo},3));
        projection_behav_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_behav_all{batchi,cond_gonogo},3));
        projection_state_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_state_all{batchi,cond_gonogo},3));
    end
end
%plot
batchi=1;fishi=9;s=1:1:length(trial_ind_session);
figure('color','w','position',[40,42,1000,1000]),kk=1;
for ii=1:4
    subplot(2,2,ii);
    for cond_gonogo=1:2; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
        end
        for cond_session=s
            x=projection_cue_all_mean{batchi,cond_gonogo}(:,cond_session);x = smoothdata(x);
            y= projection_behav_all_mean{batchi,cond_gonogo}(:,cond_session);y = smoothdata(y);
            z=projection_state_all_mean{batchi,cond_gonogo}(:,cond_session); z = smoothdata(z);
            pp=plot3( x,y,z,'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;line(cond_session)=pp;
            plot3( x(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),y(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),z(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',5,'color',clr(cond_session,:),'LineStyle',type); hold on;
            scatter3(x(1),y(1),z(1),14,clr(cond_session,:),'filled');hold on;
            n=2;scatter3(x(1:n:end),y(1:n:end),z(1:n:end),30,clr(cond_session,:),">",'filled');hold on;
        end
    end
    xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    legend([line(s)],{num2str(s')},'Location','best')
    xlabel('Cue');ylabel('Behav');zlabel('State');grid on;
    switch ii; case 1;view(125,20); case 2;view(-90,90); case 3;view(-180,0); case 4;view(-90,0);end
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
clc;clear all;
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        p=char(Path{batchi}(fishi));%'H:\1.Test US\2.Tail free——Data from 117\20210822\fish1\';
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'para.mat'));load(fullfile(p,'projection to StateCueMotion_align.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'))
        projection_state_region=[];projection_cue_region=[];projection_behav_region=[];
        for regioni=1:length(Label)
            ind_in_regioni=find(strcmp(brain_region_id(:,2),Label(regioni)));
            s=[1:length(trial_ind_session)];
            % figure, %subplot(4,2,kk);
            for cond_session=s
                for cond_gonogo=1:3 %go/nogo/all
                    avg_diff=cat(2,Difference.cue{2}(ind_in_regioni,:),Difference.behav{2}(ind_in_regioni,:),Difference.state{2}(ind_in_regioni,:));
                    w= W{2,2,2}(ind_in_regioni,:);%qr(avg_diff,0);%W{2,2,1}(ind_in_regioni,:);%
                    x=w'*squeeze(X_avg_act(ind_in_regioni,:,cond_session,cond_gonogo));
                    projection_state_region(regioni,cond_session,cond_gonogo,:)=x(3,:);
                    projection_cue_region(regioni,cond_session,cond_gonogo,:)=x(1,:);
                    projection_behav_region(regioni,cond_session,cond_gonogo,:)=x(2,:);
                    %x=W{1,2,1}(ind_in_regioni,:)'*squeeze(X(ind_in_regioni,:,cond_session,cond_gonogo));
                    %             p=plot3( x(1,:),x(2,:),x(3,:),'linewidth',1.2,'color',clr(cond_session,:),'LineStyle',type); hold on;line(cond_session)=p;
                    %             plot3( x(1,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),x(2,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),x(3,ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',2,'color',clr(cond_session,:),'LineStyle',type); hold on;
                    %             scatter3(x(1,1),x(2,1),x(3,1),14,clr(cond_session,:),'filled');hold on;
                    %             n=3;scatter3(x(1,1:n:end),x(2,1:n:end),x(3,1:n:end),14,clr(cond_session,:),">",'filled');hold on;
                end
            end
            %     legend([line(s)],[sessionx(s)],'Location','best')
            %     xlabel('Cue');ylabel('Behav');zlabel('State');grid on;title(Label(regioni));
        end
        save([p,'projection to StateCueMotion_align.mat'],'projection_state_region','projection_cue_region','projection_behav_region','-append');
    end
end

%% summary all
clc;clear all;
clr_cmap =addcolorplus(300);addcolorplus(336);%cmap = GetColormap('jet',81);
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'}; 
load(fullfile(char(Path{1}(1)),'brain_region_related_statistic.mat'),'Label')
regionxx=categorical(Label);regionxx = reordercats(regionxx,Label);
region_id=setdiff(1:length(regionxx),[18:20,23,24,27,28]);
projection_state_region_all={};projection_cue_region_all={};projection_behav_region_all={};
for batchi=1:4
    for fishi=num{batchi}
        p=char(Path{batchi}(fishi));
        load([p,'projection to StateCueMotion_align.mat'],'projection_state_region','projection_cue_region','projection_behav_region')
        projection_state_region_all{batchi,fishi}=projection_state_region;
        projection_cue_region_all{batchi,fishi}= projection_cue_region;
        projection_behav_region_all{batchi,fishi}=projection_behav_region;
    end
end
projection_cue_all=[];projection_behav_all=[];projection_state_all=[];projection_cue_all_mean={};projection_behav_all_mean={};projection_state_all_mean={};
for cond_gonogo=1%go/nogo/all
    for batchi=1:4
            kk=1;
            for fishi=num{batchi}
                NAME=projection_cue_region_all;
                aa=squeeze( NAME{batchi,fishi}(:,:,2,:));aa=preprocessfunc(aa,2);
                bb=squeeze( NAME{batchi,fishi}(:,:,1,:));bb=preprocessfunc(bb,2);
                c=(aa-repmat(aa(:,1,:),1,size(aa,2),1));%(aa-bb);abs(aa-repmat(aa(:,1,:),1,size(aa,2),1));%c=preprocessfunc(c);
                projection_cue_all(kk,:,cond_gonogo,:,:)=c;
                NAME=projection_behav_region_all;
                aa=squeeze( NAME{batchi,fishi}(:,:,1,:));aa=preprocessfunc(aa,2);
                bb=squeeze( NAME{batchi,fishi}(:,:,2,:));bb=preprocessfunc(bb,2);
                 c=(aa-repmat(aa(:,1,:),1,size(aa,2),1));%(aa-bb);c=(c-repmat(c(:,1,:),1,size(c,2),1));%c=preprocessfunc(c);
                projection_behav_all(kk,:,cond_gonogo,:,:)=c;
                NAME=projection_state_region_all;
                aa=squeeze( NAME{batchi,fishi}(:,:,3,:));aa=preprocessfunc(aa,2);
                c=abs(aa-repmat(aa(:,1,:),1,size(aa,2),1));
                projection_state_all(kk,:,cond_gonogo,:,:)=c;%squeeze(projection_state_region{typei,batchi}(:,cond_session,2,:))- squeeze(projection_state_region{typei,batchi}(:,1,2,:))+squeeze(projection_state_region{typei,batchi}(:,cond_session,1,:))-squeeze( projection_state_region{typei,batchi}(:,1,1,:));
                kk=kk+1;
            end
            projection_cue_all_mean{batchi,cond_gonogo}=squeeze(nanmean( projection_cue_all(:,:,cond_gonogo,:,:),1));
            projection_behav_all_mean{batchi,cond_gonogo}=squeeze(nanmean(projection_behav_all(:,:,cond_gonogo,:,:),1));
            projection_state_all_mean{batchi,cond_gonogo}=squeeze(nanmean(projection_state_all(:,:,cond_gonogo,:,:),1));
    end
end
%% plot cue    
figure('position',[1 1 1800 1000]); tiledlayout(4,8,'TileSpacing','Compact');
for batchi=1:4;
    for cond_gonogo=1
        b=projection_cue_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,1);
        a=[-2 10];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(hot);colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
        end
    end
end
%plot behav
figure('position',[1 1 1800 1000]); tiledlayout(4,8,'TileSpacing','Compact');
for batchi=1:4;
    for cond_gonogo=1
        b=projection_behav_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,1);
        a=[-0.1 0.6];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(hot);colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
        end
    end
end
%plot state
figure('position',[1 1 1800 1000]); tiledlayout(4,8,'TileSpacing','Compact');
for batchi=1:4;
    for cond_gonogo=1
        b=projection_state_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,1);
        a=[-0 10];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(hot);colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
        end
    end
end
%% region_trajectory
load(fullfile(char(Path{1}(1)),'para.mat'),'frame');
load([p,'projection to StateCueMotion_align.mat'],'bin');
for batchi=1:4
    figure('color','w','position',[1 1 800 1000]),tiledlayout(3,2,'TileSpacing','Compact');
    for regioni=[3,5,7,25,26,29];
        for cond_gonogo=1; %go/nogo/all
            switch cond_gonogo
                case 1
                    type='-';
                case 2
                    type='--';
            end
            nexttile;
            s=1:size(projection_cue_all_mean{batchi,cond_gonogo},2);
            clr=hsv(length(s));
            for cond_session=s
                x=squeeze(projection_cue_all_mean{batchi,cond_gonogo}(regioni,:,:))';x=x(:,cond_session);x = smoothdata(x);%x=preprocessfunc(x);
                y= squeeze(projection_behav_all_mean{batchi,cond_gonogo}(regioni,:,:))';y=y(:,cond_session);y = smoothdata(y);%y=preprocessfunc(y);
                z=squeeze(projection_state_all_mean{batchi,cond_gonogo}(regioni,:,:))';z=z(:,cond_session);z = smoothdata(z);%z=preprocessfunc(z);
                
                pp=plot3( x,y,z,'linewidth',2,'color',clr(cond_session,:),'LineStyle', type); hold on;line(cond_session)=pp;
                plot3( x(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),y(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),z(ceil(frame.cs_start/bin):ceil(frame.cs_end/bin)),'linewidth',5,'color',clr(cond_session,:),'LineStyle', type); hold on;
                scatter3(x(1),y(1),z(1),14,clr(cond_session,:),'filled');hold on;
                n=2;scatter3(x(1:n:end),y(1:n:end),z(1:n:end),30,clr(cond_session,:),">",'filled');hold on;
            end
        end
       % legend([line(s)],{num2str(s')},'Location','best');
        title(Label{regioni});view(0,90)
        xlabel('Cue');ylabel('Behav');zlabel('State');grid on;
    end
end
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