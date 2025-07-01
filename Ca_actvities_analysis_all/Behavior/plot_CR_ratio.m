mat_name=[];
path={['H:\1.Test US\2.Tail free！！Data from 117\20220709\fish1\',mat_name],...%no Cb after Learning
    ['H:\1.Test US\2.Tail free！！Data from 117\20220802\fish2\',mat_name],...%not typical
    ['H:\1.Test US\2.Tail free！！Data from 117\20220803\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220814\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220915\fish1\',mat_name],...%calcium wrong
    ['H:\1.Test US\2.Tail free！！Data from 117\20210402\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210709\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210805\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20230309\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20230602\fish1\',mat_name],...  %音頁learner��
    ['H:\1.Test US\2.Tail free！！Data from 117\20230807\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210514\fish2\',mat_name],...%un shuffle
    ['H:\1.Test US\2.Tail free！！Data from 117\20210702\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220831\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20221108\fish1\',mat_name ],...% wrong
    ['H:\1.Test US\2.Tail free！！Data from 117\20230309\fish2\',mat_name ] %1
    };%learner
for batchi=1:length(path)
if isfolder(fullfile(path{batchi},'behavior\'))
    path{batchi}
    [status,msg] = movefile(fullfile(path{batchi},'behavior\'),fullfile(path{batchi},'behav\'))
    status
end
end
load(fullfile(path{1},'para.mat'));
trial_ind_session=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            trial_ind_session(ii,:)=[trial.hab(2),trial.hab(3)];
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            trial_ind_session(ii,:)=[(trial.hab(3)+(ii-2)*trial.acq_block_trial)+1,(trial.hab(3)+(ii-1)*trial.acq_block_trial)];
        case trial.acq_block_num+2
            trial_ind_session(ii,:)=[trial.test(2),trial.test(3)];
    end
end
trial_ind_session_slid=[];
slidwin=6;
trial_ind_session_slid(1,:)=[trial.hab(2),trial.hab(3)];kk=2;
for ii=trial.acq(2):trial.acq(3)-slidwin+1
    trial_ind_session_slid(kk,:)=[ii, ii+slidwin-1];kk=kk+1;
end
trial_ind_session_slid(kk,:)=[trial.test(2),trial.test(3)];

CR_ratio_1=[];CR_ratio_2=[];CR_ratio_fishname=string;
for jj=1:length(path)
    p=path{jj};
    load(fullfile(p,'para.mat'));
    load(fullfile(p,'behav\behavior_from_Results_of_alltheta.mat'));
    nn=[p(end-14:end-7),p(end-5:end-1)];
    CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.2*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
    for ii=1:size(trial_ind_session,1)
        CR_ratio_1(ii,jj)=length(unique(CR(find(CR(:,1)<=trial_ind_session(ii,end) & CR(:,1)>=trial_ind_session(ii,1)),1)))/( trial_ind_session(ii,end)-trial_ind_session(ii,1)+1);
    end
    for ii=1:size(trial_ind_session_slid,1)
        CR_ratio_2(ii,jj)=length(unique(CR(find(CR(:,1)<=trial_ind_session_slid(ii,end) & CR(:,1)>=trial_ind_session_slid(ii,1)),1)))/(trial_ind_session_slid(ii,end)-trial_ind_session_slid(ii,1)+1);
    end
    CR_ratio_fishname(jj,1)=nn;
end
figure,plot(mean(CR_ratio_2,2))
figure,plot(mean(CR_ratio_1,2))