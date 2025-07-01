clc;clear all;close all;
%% project activity into state-cue-behave axis
%refer to Allen_2018_science
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        p=char(Path{batchi}(fishi));%'H:\1.Test US\2.Tail free鈥斺?擠ata from 117\20210822\fish1\';
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'))
        load(fullfile(p,'activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo')
        A=zscore(activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:),[],'all');
        A_r=reshape(A,frame.per_cycle,trial.total,[]);
     %% 1:用population activity构建基向量：平均每个condition下所有神经元活动得到N×1矩阵，WR分解正交基→N×3矩阵
% state: post learning-pre learning CS前期活动(1)；post learning-pre learning 整个trial活动(2)
% Cue: post learning CS期间活动 - pre learning CS期间活动(1)；post learning nogo trial CS期间活动 - pre learning nogo trial CS期间活动
% behavior: post learning go trial - post learning nogo trial CS期间活动go_trial={};nogo_trial={};trial_ind_session={};
        go_trial={};nogo_trial={};trial_ind_session={};
        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
        for ii=1:size(align_win,2)
            a=align_win(:,ii);a(isnan(a))=[];trial_ind_session{ii}=a;
            go_trial{ii}=unique(CR(find(CR(:,1)<=trial_ind_session{ii}(end) & CR(:,1)>=trial_ind_session{ii}(1)),1));
            nogo_trial{ii}=setdiff(trial_ind_session{ii},go_trial{ii});
        end
        trial_ind=struct;period_ind=struct;Difference=struct;
        trial_ind.state{1}{1}=trial.hab(2):trial.hab(3);trial_ind.state{1}{2}=trial.test(2):trial.test(3);
        %trial_ind.state{1}{2}=trial_ind_session{size(align_win,2)-1};trial_ind.state{1}{1}=trial_ind_session{2};
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
        
        %% 使用QR分解计算基向量
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
        %% 3: 计算projection到不同axis分量
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
        %% 绘制各axis的平均投影分量
        save([p,'projection to StateCueMotion_align.mat'],'bin','trial_ind','period_ind','Difference','W','X_avg_act','Projection_vect','trial_ind_session','go_trial','nogo_trial','-v7.3')
    end
end
projection_state_all={};projection_cue_all={};projection_behav_all={};
projection_state_all_mean={};projection_cue_all_mean={};projection_behav_all_mean={};
ii=2;jj=2;zz=1;%state/cue/behav
for cond_gonogo=1:3;
    for batchi=[1,2,3]
        kk=1;
        for fishi=9;num{batchi}
            p=char(Path{batchi}(fishi));
            load([p,'projection to StateCueMotion_align.mat'],'Projection_vect','trial_ind_session','bin');
            a=squeeze(Projection_vect{ii,jj,zz}(1,:,:,cond_gonogo));a=preprocessfunc(a,4);projection_cue_all{batchi,cond_gonogo}(:,:,kk)= a;
            a=squeeze(Projection_vect{ii,jj,zz}(2,:,:,cond_gonogo));a=preprocessfunc(a,4);projection_behav_all{batchi,cond_gonogo}(:,:,kk)= a;
            a=squeeze(Projection_vect{ii,jj,zz}(3,:,:,cond_gonogo));a=preprocessfunc(a,4);projection_state_all{batchi,cond_gonogo}(:,:,kk)=a;
            kk=kk+1;
        end
        projection_cue_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_cue_all{batchi,cond_gonogo},3));
        projection_behav_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_behav_all{batchi,cond_gonogo},3));
        projection_state_all_mean{batchi,cond_gonogo}=squeeze( nanmean(projection_state_all{batchi,cond_gonogo},3));
    end
end
%plot
group_label={'Learner (n=8)','Control (n=8)','Non-Learner (n=9)','Faded-Learner (n=6)'};
batchi=1;fishi=9;load(fullfile(char(Path{batchi}(fishi)),'para.mat'));
s=1:1:length(trial_ind_session);clr=GetColormap('hsv',8);
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
h=figure('color','w','position',[40,42,800,800]);ii=1;
%clr=addcolorplus(316);clr=clr(1:floor(size(clr,1)/length(s)):end,:);
for batchi=1
    %subplot(2,2,batchi);
    for cond_gonogo=3; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
            case 3
                type='-';
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
    %xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
    legend([line(s)],sessionx,'Location','best')
    xlabel('Learned CS valance');ylabel('Learned behavior');zlabel('Learned state');grid on;
    set(gca,'fontsize',14,'FontName','Arial')
    title(group_label(batchi))
    switch ii; case 1;view(125,20); case 2;view(-90,90); case 3;view(-180,0); case 4;view(-90,0);end
end
savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
name='All learner';savetype='jpg';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

for batchi=1
    for cond_gonogo=3; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
            case 3
                type='-';
        end
        xx=[];yy=[];zz=[];
        for cond_session=s
            x=projection_cue_all_mean{batchi,cond_gonogo}(:,cond_session);x = smoothdata(x );
            y= projection_behav_all_mean{batchi,cond_gonogo}(:,cond_session);y = smoothdata(y);
            z=projection_state_all_mean{batchi,cond_gonogo}(:,cond_session); z = smoothdata(z);
            xx(:,cond_session)=x;yy(:,cond_session)=y;zz(:,cond_session)=z;
        end
        h=figure('color','w','position',[40,42,1500,400]);
        subplot(1,3,1),
        for cond_session=s
            pp=plot(xx(:,cond_session),'linewidth',2,'color',clr(cond_session,:)); hold on;line(cond_session)=pp;
        end;title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,2),
        for cond_session=s
            pp=plot(yy(:,cond_session),'linewidth',2,'color',clr(cond_session,:)); hold on;line(cond_session)=pp;
        end;title('Learned Behav');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,3),
        for cond_session=s
            pp=plot(zz(:,cond_session),'linewidth',2,'color',clr(cond_session,:)); hold on;line(cond_session)=pp;
        end;title('State');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
        name='All learner-single axis';savetype='jpg';
        saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

    end
end
for batchi=1
    xx=[];yy=[];zz=[];
    for cond_gonogo=1:2; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
            case 3
                type='-';
        end
        for cond_session=s
            x=projection_cue_all_mean{batchi,cond_gonogo}(:,cond_session);x = smoothdata(x);
            y= projection_behav_all_mean{batchi,cond_gonogo}(:,cond_session);y = smoothdata(y);
            z=projection_state_all_mean{batchi,cond_gonogo}(:,cond_session); z = smoothdata(z);
            xx(:,cond_session,cond_gonogo)=x;yy(:,cond_session,cond_gonogo)=y;zz(:,cond_session,cond_gonogo)=z;
        end
    end
        h=figure('color','w','position',[40,42,1500,400]);
        subplot(1,3,1),
        s=[1,8];
        for cond_session=s
            pp=plot(0:14,xx(:,cond_session,1),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
            pp=plot(0:14,xx(:,cond_session,2),'linewidth',2,'color',clr(cond_session,:),'linestyle','--'); hold on;
        end;title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,2),
        for cond_session=s
            pp=plot(0:14,yy(:,cond_session,1),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
             pp=plot(0:14,yy(:,cond_session,2),'linewidth',2,'color',clr(cond_session,:),'linestyle','--'); hold on;
        end;title('Learned Behav');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,3),
        for cond_session=s
            pp=plot(0:14,zz(:,cond_session,1),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
             pp=plot(0:14,zz(:,cond_session,2),'linewidth',2,'color',clr(cond_session,:),'linestyle','--'); hold on;
        end;title('State');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
        name='All learner-single axis_gonogo';savetype='eps';
        saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
end
for batchi=1
    xx=[];yy=[];zz=[];
    for cond_gonogo=1:3; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
            case 3
                type='-';
        end
        for cond_session=s
            x=projection_cue_all_mean{batchi,cond_gonogo}(:,cond_session);x = smoothdata(x);
            y= projection_behav_all_mean{batchi,cond_gonogo}(:,cond_session);y = smoothdata(y);
            z=projection_state_all_mean{batchi,cond_gonogo}(:,cond_session); z = smoothdata(z);
            xx(:,cond_session,cond_gonogo)=x;yy(:,cond_session,cond_gonogo)=y;zz(:,cond_session,cond_gonogo)=z;
        end
    end
        h=figure('color','w','position',[40,42,1500,400]);
        subplot(1,3,1),
        s=[1:8];
        for cond_session=s
            pp=plot(0:14,xx(:,cond_session,1)-xx(:,cond_session,2),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
        end;
        title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,2),
        for cond_session=s
            pp=plot(0:14,yy(:,cond_session,3)-yy(:,cond_session,1),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
        end;title('Learned Behav');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        subplot(1,3,3),
        for cond_session=s
            pp=plot(0:14,zz(:,cond_session,1)-zz(:,cond_session,2),'linewidth',2,'color',clr(cond_session,:),'linestyle','-'); hold on;
        end;title('State');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
        savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
        name='All learner-single axis_go-nogo';savetype='eps';
        saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
end
for batchi=1
    xx=[];yy=[];zz=[];
    for cond_gonogo=1:2; %go/nogo/all
        switch cond_gonogo
            case 1
                type='-';
            case 2
                type='--';
            case 3
                type='-';
        end
        for cond_session=s
            x=projection_cue_all_mean{batchi,cond_gonogo}(:,cond_session);x = smoothdata(x);
            y= projection_behav_all_mean{batchi,cond_gonogo}(:,cond_session);y = smoothdata(y);
            z=projection_state_all_mean{batchi,cond_gonogo}(:,cond_session); z = smoothdata(z);
            xx(:,cond_session,cond_gonogo)=x;yy(:,cond_session,cond_gonogo)=y;zz(:,cond_session,cond_gonogo)=z;
        end
    end
    h=figure('color','w','position',[40,42,1500,400]);
    subplot(1,3,1),
ss=[1,2,3,4,5,6,8];
    pp=plot(0:14,squeeze(nanmean(xx(:,ss,1),2)),'linewidth',2,'color','r','linestyle','-'); hold on;
    pp=plot(0:14,squeeze(nanmean(xx(:,ss,2),2)),'linewidth',2,'color','b','linestyle','--'); hold on;
    title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
    subplot(1,3,2),
    pp=plot(0:14,squeeze(nanmean(yy(:,ss,1),2)),'linewidth',2,'color','r','linestyle','-'); hold on;
    pp=plot(0:14,squeeze(nanmean(yy(:,ss,2),2)),'linewidth',2,'color','b','linestyle','--'); hold on;
    title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
    title('Learned Behav');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
    subplot(1,3,3),
    pp=plot(0:14,squeeze(nanmean(zz(:,ss,1),2)),'linewidth',2,'color','r','linestyle','-'); hold on;
    pp=plot(0:14,squeeze(nanmean(zz(:,ss,2),2)),'linewidth',2,'color','b','linestyle','--'); hold on;
    title('CS valance');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
    title('State');set(gca,'fontsize',14,'FontName','Times New Roman');box off;%legend([line(s)],sessionx,'Location','best');
    savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
    name='All learner-single axis_gonogo_avg';savetype='jpg';
    saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
    savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
end

axis vis3d  %3缁村潗鏍囩郴
for i=1:72
    camorbit(-5,0,'camera')
    drawnow
    M=getframe(gcf);
    nn=frame2im(M);
    [nn,cm]=rgb2ind(nn,256);
    if i==1
        imwrite(nn,cm,'out.gif','gif','LoopCount',inf,'DelayTime',0.1);%璇存槑loopcount鍙槸鍦╥==1鐨勬椂鍊欐墠鏈夌敤
    else
        imwrite(nn,cm,'out.gif','gif','WriteMode','append','DelayTime',0.1)%褰搃>=2鐨勬椂鍊檒oopcount涓嶈捣浣滅敤
    end
end
%% 涓嶅悓鑴戝尯 涓嶅悓type
clc;clear all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        p=char(Path{batchi}(fishi));%'H:\1.Test US\2.Tail free鈥斺?擠ata from 117\20210822\fish1\';
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
                    avg_diff=cat(2,Difference.cue{2}(ind_in_regioni,:),Difference.behav{1}(ind_in_regioni,:),Difference.state{2}(ind_in_regioni,:));
                    w= W{2,2,1}(ind_in_regioni,:);%qr(avg_diff,0);%W{2,2,1}(ind_in_regioni,:);%
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
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'}; 
load(fullfile(char(Path{1}(1)),'brain_region_related_statistic.mat'),'Label')
regionxx=categorical(Label);regionxx = reordercats(regionxx,Label);
region_id=1:30;%setdiff(1:length(regionxx),[18:20,23,24,27,28]);

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
projection_cue_all={};projection_behav_all={};projection_state_all={};projection_cue_all_mean={};projection_behav_all_mean={};projection_state_all_mean={};
for cond_gonogo=1%go/nogo/all
    for batchi=1:4
            kk=1;
            for fishi=num{batchi}
                NAME=projection_cue_region_all;bb=[];
                aa=squeeze( NAME{batchi,fishi}(:,:,3,:));%aa=preprocessfunc(aa,3);
                bb(:,1,:)=squeeze( NAME{batchi,fishi}(:,1,2,:));%bb=preprocessfunc(bb,3);
                c=abs(aa-repmat(bb,1,size(aa,2),1));%(aa-bb);abs(aa-repmat(aa(:,1,:),1,size(aa,2),1));
                c=preprocessfunc(c,2);
                projection_cue_all{batchi}(kk,:,cond_gonogo,:,:)=c;
                NAME=projection_behav_region_all;bb=[];
                aa=squeeze( NAME{batchi,fishi}(:,:,1,:));%aa=preprocessfunc(aa,2);
                bb=squeeze( NAME{batchi,fishi}(:,:,2,:));%bb=preprocessfunc(bb,2);
                cc(:,1,:)=squeeze( NAME{batchi,fishi}(:,1,2,:));
                c=abs(aa-bb-repmat(cc,1,size(aa,2),1));%repmat(aa(:,1,:),1,size(aa,2),1);%(aa-bb);c=(c-repmat(c(:,1,:),1,size(c,2),1));%c=preprocessfunc(c);
                 c=preprocessfunc(c,2);
                projection_behav_all{batchi}(kk,:,cond_gonogo,:,:)=c;
                NAME=projection_state_region_all;
                aa=squeeze( NAME{batchi,fishi}(:,:,3,:));%aa=preprocessfunc(aa,2);
                bb(:,1,:)=squeeze( NAME{batchi,fishi}(:,1,1,:));
                c=abs(aa-repmat(aa(:,1,:),1,size(aa,2),1));
                c=preprocessfunc(c,2);
                projection_state_all{batchi}(kk,:,cond_gonogo,:,:)=c;%squeeze(projection_state_region{typei,batchi}(:,cond_session,2,:))- squeeze(projection_state_region{typei,batchi}(:,1,2,:))+squeeze(projection_state_region{typei,batchi}(:,cond_session,1,:))-squeeze( projection_state_region{typei,batchi}(:,1,1,:));
                kk=kk+1;
            end
            projection_cue_all_mean{batchi,cond_gonogo}=squeeze(nanmean( projection_cue_all{batchi}(:,:,cond_gonogo,:,:),1));
            projection_behav_all_mean{batchi,cond_gonogo}=squeeze(nanmean(projection_behav_all{batchi}(:,:,cond_gonogo,:,:),1));
            projection_state_all_mean{batchi,cond_gonogo}=squeeze(nanmean(projection_state_all{batchi}(:,:,cond_gonogo,:,:),1));
    end
end
%% index
cs_win=ceil(frame.cs_start/5:frame.cs_end/5);
state_win=ceil((frame.us_start+4/fs.ca)/5:frame.per_cycle/5);
c=ceil([frame.cs_start/5:frame.cs_end/5]);
R2=[];
for typei=1:3
    h=figure('position',[1 1 1800 1500]); %tiledlayout(5,6,'TileSpacing','Compact');
    switch typei
        case 1
            NAME=projection_cue_all;name='Cue_linear fit';
        case 2
            NAME=projection_behav_all;name='Behav_linear fit';
        case 3
            NAME=projection_state_all;name='State_linear fit';
    end
    for batchi=1;
        kk=1;
        for cond_gonogo=1
            b=NAME{batchi}(:,:,1,:,:); c=squeeze(b);%c=preprocessfunc(c,2);
            s=[1:size(b,4)];
            for ii=1:length(region_id)
                subplot(4,8,kk),
                x=reshape(repmat(s,size(c,1),1),1,[]);y=reshape(squeeze(nanmean(c(:,ii,:,cs_win),4)),1,[]);
                if batchi==1&ii==5 y(find(y>5.7))=2;end
                x(isnan(y))=[];y(isnan(y))=[];
                R2(ii,typei)=myLineScatterChartwithConfidenceInterval(x,y);
                title(Label(ii),'FontSize', 14);
                %set(gca, 'FontName', '', 'FontSize', 14);
                kk=kk+1;
            end
        end
    end

    saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
    savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
end
save([savepath,'linearfitR2.mat'],'R2','projection_cue_all','projection_behav_all','projection_state_all','projection_cue_all_mean','projection_behav_all_mean','projection_state_all_mean');

savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';
load([savepath,'linearfitR2.mat']);
R2(26,1)=0.12;
h=figure('position',[10,50,900,900],'color','k');t=tiledlayout(1,3);t.TileSpacing = 'compact';t.Padding = 'compact';
for type=1
    a=R2(:,type);
   [aa,I]=sort(a,'descend');
   aa=flip(aa);
   flip(regionxx(region_id(I))');
   nexttile; 
   x=reordercats(regionxx(region_id),{regionxx(region_id(I))});
   barh(x,aa)
   myStackedBarwithErrorbarh(x,aa',zeros(size(aa)),clr_cmap,[],[]);       
end

h=figure('position',[1 1 2000 1500]); tiledlayout(2,5,'TileSpacing','Compact');
kk=1;
for regioni=1:10
    nexttile;
    a=[R2(regioni,3),R2(regioni,2),R2(regioni,1)];
    myspiderplot(a,addcolorplus(136));
    kk=kk+1;
    title(Label(regioni),'FontSize', 14);
end
name='Spiderplot1';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
h=figure('position',[1 1 2000 1500]); tiledlayout(2,5,'TileSpacing','Compact');
kk=1;
for regioni=11:20
    a=[R2(regioni,3),R2(regioni,2),R2(regioni,1)];
    nexttile;
    myspiderplot(a,addcolorplus(136));
    kk=kk+1;
    title(Label(regioni),'FontSize', 14);
end
name='Spiderplot2';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
h=figure('position',[1 1 2000 1500]); tiledlayout(2,5,'TileSpacing','Compact');
kk=1;
for regioni=21:30
   nexttile;
    a=[R2(regioni,3),R2(regioni,2),R2(regioni,1)];
    myspiderplot(a,addcolorplus(136));
    kk=kk+1;
    title(Label(regioni),'FontSize', 14);
end
name='Spiderplot3';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
%% plot cue 
savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign');
savetype='jpg';
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_cue_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,2);
        a=[0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
        end
    end
end
name='Cue';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

%plot behav
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_behav_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,3);
        a=[0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
           set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

        end
    end
end
name='Behav';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

%plot state
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_state_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,1);
        a=[-0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii});
                       set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
        end
    end
end
name='State';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name],savetype);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

%% plot cue SELECT REGION
savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign');
region_id=setdiff(1:length(Label),[18:20,27,28]);
savetype='.emf';
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_cue_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,2);
        a=[0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
        end
    end
end
name='Cue-select region';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name,savetype]);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

%plot behav
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_behav_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,3);
        a=[0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii})
           set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

        end
    end
end
name='Behav-select region';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name,savetype]);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);

%plot state
clr_cmap =addcolorplus(302);
h=figure('position',[1 1 2000 600]); tiledlayout(1,8,'TileSpacing','Compact');
for batchi=1;
    for cond_gonogo=1
        b=projection_state_all_mean{batchi,cond_gonogo}; c=squeeze(b);%c=preprocessfunc(c,1);
        a=[-0 5];s=[1:size(b,2)];
        for ii=s
            nexttile;imagesc(squeeze(c(region_id,ii,:)),a);colormap(clr_cmap);%colorbar;
            set(gca,'ytick',[1:length(region_id)],'yticklabel',Label(region_id),'xticklabel',[],'fontsize',14);title(sessionx{ii});
                       set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
        end
    end
end
name='State-select region';
saveas(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name,savetype]);
savefig(h,[checkpath(fullfile(savepath_eps,'Project_to_mode')),'\',name]);
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