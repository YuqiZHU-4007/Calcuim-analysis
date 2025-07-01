clc;clear all;
set(0,'defaultfigurecolor','w');

mat_name=[];
batchpath={['H:\1.Test US\2.Tail free！！Data from 117\20220709\fish1\',mat_name],...%no Cb after Learning 
    ['H:\1.Test US\2.Tail free！！Data from 117\20220802\fish2\',mat_name],...%not typical
    ['H:\1.Test US\2.Tail free！！Data from 117\20220803\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220814\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220915\fish1\',mat_name],...%calcium wrong
    ['H:\1.Test US\2.Tail free！！Data from 117\20210402\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210709\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210805\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20230309\fish1\',mat_name],...     %['H:\1.Test US\2.Tail free！！Data from 117\20230602\fish1\',mat_name],...  %音頁learner��
    ['H:\1.Test US\2.Tail free！！Data from 117\20230807\fish1\',mat_name]
    };%learner

bin=32;
for batchi=1:length(batchpath)
    load(fullfile(batchpath{batchi},'para.mat'));
    a=load(fullfile(batchpath{batchi},'behav\behavior_from_Results_of_alltheta.mat'));
    sliding_win=reshape([1:ceil(frameb.per_cycle/bin)*bin],bin,[]);
    CS_bin=[frameb.cs_start/bin frameb.cs_end/bin];
    US_bin=frameb.us_start/bin;
    freq=zeros(size(sliding_win,2),trial.total);
    for ii=1:trial.total
        ind=a.re_startpoint_sd(find(a.re_startpoint_sd(:,1)==ii),2);
        if ~isempty(ind)
            [C,~,ib]=intersect(ind,sliding_win);
            [I,J] = ind2sub(size(sliding_win),ib);%sliding_win(I,J)
            for jj=unique(I)'
                ind_jj=find(I==jj);
                freq(J(ind_jj),ii)=length(ind_jj);%/(bin/fs.ca);
            end
        end
    end
    figure,
    for jj=1:3
        ax=subplot(3,1,jj);
        switch jj
            case 1
                trial_ind=trial.hab(2):trial.hab(3); l='Hab.';
            case 2
                trial_ind=trial.acq(2):trial.acq(3); l='Acq.';
            case 3
                trial_ind=trial.test(2):trial.test(3);l='Tst.';
        end
        p1=plot(freq(:,trial_ind),'color',[0.5 0.5 0.5]);hold on;        
        %CS bar
        color=[1 0 0];
        patch1=patch([CS_bin(1),CS_bin(2),CS_bin(2),CS_bin(1)],...
            [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],...
            color,'edgecolor',color,'facecolor',color,'edgealpha',0.3,'facealpha',0.3);hold on %CS on-off
        l1=line([US_bin US_bin],[ax.YLim(1) ax.YLim(2)],'color',[0 0 0],'linestyle','--','linewidth',1.2);hold on
        
        p1=plot(freq(:,trial_ind),'color',[0.5 0.5 0.5]);hold on;
        p2=plot(mean(freq(:,trial_ind),2),'r','linewidth',2);box on
        legend([p1(1),p2,patch1,l1],{'Each trial','Averaged','CS','US'});%legend('boxoff')
        ylabel('Number of shaking','fontsize',12);title(l,'fontsize',12);
    end
end