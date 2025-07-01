%1:trial num;2:start trial;3:end trial
trial.spon_bef=[6*3 1 6*3];
trial.hab = [15 trial.spon_bef(3)+1 trial.spon_bef(3)+15];
trial.acq_block_num=4;
trial.acq_block_trial=6;
trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
trial.test =[6 trial.acq(3)+1 trial.acq(3)+6];
trial.spon_aft=[6*3 trial.test(3)+1 trial.test(3)+6*3];
trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);
 
a=re_startpoint_sd;
a(find(re_startpoint_sd(:,1)<trial.hab(2)),:)=0;
a(find(re_startpoint_sd(:,1)>trial.hab(3)),:)=0;
a(find(re_startpoint_sd(:,2)<600),:)=0;
a(find(re_startpoint_sd(:,2)>60*14.8),:)=0;
a(find(a(:,1)==0),:)=[];
a(:,1)=a(:,1)-trial.hab(2)+1;