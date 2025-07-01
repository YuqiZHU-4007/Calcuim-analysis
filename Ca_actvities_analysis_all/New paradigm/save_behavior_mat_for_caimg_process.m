clc;clear all
global fs
global frameb
global trial
[fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=10;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num)
T_non_spon=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frameb.per_cycle+1:(trial.hab(3))*frameb.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(2)-1)*frameb.per_cycle+1:(trial.test(3))*frameb.per_cycle;
    end
     T_non_spon=[T_non_spon,t];
end
T_spon={};
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.spon_bef(2)-1)*frameb.per_cycle+1:(trial.spon_bef(3))*frameb.per_cycle;
        case mat2cell([2:trial.acq_block_num],1,ones(1,trial.acq_block_num-1))
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
        case trial.acq_block_num+1
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(3))*frameb.per_cycle+1:(trial.total)*frameb.per_cycle;
    end
    T_spon{ii}=t;
end

is_addiomission=0;
%%
pp=uigetdir('H:\1.Test US\2.Tail free——Data from 117\20210805\fish1\behav');
load(fullfile(pp,'Results_of_alltheta.mat'))
% if is_addiomission
%     alltheta_omission=198001:199800;
%     alltheta=alltheta(setdiff([1:length(alltheta)],alltheta_omission));
% end
behav=alltheta(T_non_spon);behav_spon=[];
for ii=1:trial.acq_block_num+2
behav_spon(:,ii)=alltheta(T_spon{ii});
end
save([pp '\behav'],'behav');save([pp '\behav_spon.mat'],'behav_spon');

re_startpoint_sd=[];
y_3sd=alltheta(T_non_spon);
startpoint_sd=env.env_loc;
re_startpoint_sd(:,1)=fix(startpoint_sd/frameb.per_cycle)+1;
re_startpoint_sd(:,2)=mod(startpoint_sd,frameb.per_cycle);
save([pp '\behavior_from_Results_of_alltheta'],'y_3sd','startpoint_sd','re_startpoint_sd');
a=0;trial.spon_bef=[a min(1,a) a];
a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
trial.acq_block_num=8;
trial.acq_block_trial=3;
trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);
[a,fishi]=fileparts(fileparts(path));[~,datei]=fileparts(a);
[p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint_sd,fishi,datei,path,1);


% a=0;trial.spon_bef=[a min(1,a) a];
% a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
% trial.acq_block_num=8;
% trial.acq_block_trial=3;
% trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
% a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
% a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
% trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);
% 
% pp='H:\2.Summary of behavior\行为统计\';
% list1=dir(pp); CR_ratio={};CR_onset={};
% for ll1=3:length(list1)
%     pp1= (fullfile(list1(ll1).folder,list1(ll1).name));
%     list2=dir(pp1);kk=1;
%     for ll2=3:length(list2)
%         pp2=fullfile(list2(ll2).folder,list2(ll2).name);
%         path=list2(ll2).folder;
%         name=list2(ll2).name;
%         if isfolder(pp2)
%             load(fullfile(pp2,'behavior_from_Results_of_alltheta.mat'))
%             ind=zeros(1,10);ind_fraction=zeros(1,10);onset=nan(1,10);
%             for jj=1:10
%                 switch jj
%                     case 1
%                         trial_ind=trial.hab(2):trial.hab(3);
%                         l{jj}='Pre cond.';
%                     case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
%                         trial_ind=trial.hab(3)+(jj-2)*trial.acq_block_trial+1 :(trial.hab(3)+(jj-1)*trial.acq_block_trial);
%                         l{jj}=['Cond.' num2str(jj-1)];
%                     case 10
%                         trial_ind=trial.test(2):trial.test(3);
%                         l{jj}='Post cond.';
%                 end
%                 b=(find( re_startpoint_sd(:,1)>=trial_ind(1) & re_startpoint_sd(:,1)<=trial_ind(end) & re_startpoint_sd(:,2)>=frameb.cs_start &  re_startpoint_sd(:,2)<=frameb.us_start));
%                 a=re_startpoint_sd(b,:);
%                 ind(jj)=length(unique(a(:,1)));
%                 ind_fraction(jj)=ind(jj)./length(trial_ind);oo=[];kkk=1;
%                 for zz=unique(a(:,1))'
%                     oo(kkk,1)=min(a(find(a(:,1)==zz),2));kkk=kkk+1;
%                 end
%                 onset(jj)=mean(oo);
%             end
%             CR_ratio{ll1}(:,kk)=ind_fraction;
%             CR_onset{ll1}(:,kk)=onset;
%             kk=kk+1;
%         end
%     end
% end
% 
% a=CR_onset{4}; a=(a-frameb.cs_start)./fs.behavior;
% a=flip(a);
% 
% bin=32;
% list1=dir(pp); CR_ratio={};CR_onset={};
% for ll1=3:length(list1)
%     pp1= (fullfile(list1(ll1).folder,list1(ll1).name));
%     list2=dir(pp1);kk=1;
%     for ll2=3:length(list2)
%         pp2=fullfile(list2(ll2).folder,list2(ll2).name);
%         path=list2(ll2).folder;
%         name=list2(ll2).name;
%         if isfolder(pp2)
%             a=load(fullfile(pp2,'behavior_from_Results_of_alltheta.mat'));
%             sliding_win=reshape([1:ceil(frameb.per_cycle/bin)*bin],bin,[]);
%             CS_bin=[frameb.cs_start/bin frameb.cs_end/bin];
%             US_bin=frameb.us_start/bin;
%             freq=zeros(size(sliding_win,2),trial.total);
%             for ii=1:trial.total
%                 ind=a.re_startpoint_sd(find(a.re_startpoint_sd(:,1)==ii),2);
%                 if ~isempty(ind)
%                     [C,~,ib]=intersect(ind,sliding_win);
%                     [I,J] = ind2sub(size(sliding_win),ib);%sliding_win(I,J)
%                     for jj=unique(I)'
%                         ind_jj=find(I==jj);
%                         freq(J(ind_jj),ii)=length(ind_jj);%/(bin/fs.ca);
%                     end
%                 end
%             end
%             figure,
%             for jj=1:3
%                 ax=subplot(3,1,jj);
%                 switch jj
%                     case 1
%                         trial_ind=trial.hab(2):trial.hab(3); l='Hab.';
%                     case 2
%                         trial_ind=trial.acq(2):trial.acq(3); l='Acq.';
%                     case 3
%                         trial_ind=trial.test(2):trial.test(3);l='Tst.';
%                 end
%                 p1=plot(freq(:,trial_ind),'color',[0.5 0.5 0.5]);hold on;
%                 %CS bar
%                 color=[1 0 0];
%                 patch1=patch([CS_bin(1),CS_bin(2),CS_bin(2),CS_bin(1)],...
%                     [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],...
%                     color,'edgecolor',color,'facecolor',color,'edgealpha',0.3,'facealpha',0.3);hold on %CS on-off
%                 l1=line([US_bin US_bin],[ax.YLim(1) ax.YLim(2)],'color',[0 0 0],'linestyle','--','linewidth',1.2);hold on
%                 
%                 p1=plot(freq(:,trial_ind),'color',[0.5 0.5 0.5]);hold on;
%                 p2=plot(mean(freq(:,trial_ind),2),'r','linewidth',2);box on
%                 legend([p1(1),p2,patch1,l1],{'Each trial','Averaged','CS','US'});%legend('boxoff')
%                 ylabel('Number of shaking','fontsize',12);title(l,'fontsize',12);
%             end
%         end
%     end
% end

