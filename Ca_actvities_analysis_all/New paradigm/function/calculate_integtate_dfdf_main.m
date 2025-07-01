%计算所有roi的 integrate dfdf
%20190327-zyq
%include <calculate_integrate_dfdf.m>
function [rawacti,area]=calculate_integtate_dfdf_main(activities,cut_move,trial,frame,fs,ind_cut_trial)
%计算曲线下积分
if trial.hab(1)<=trial.test(1)
    indcut1=0;
    indcut2=trial.test(1)-trial.hab(1);
else
    indcut1=trial.hab(1)-trial.test(1);
    indcut2=0;
end
area_win_acq=frame.cs_start:frame.us_start-1;%frame.us_start
area_win_hab_tst=frame.cs_start:frame.cs_end-1;
type_cal_area='event';%type={'lowest','lowest_all','polyfit','first','mean','sum','event'};
ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;%取取时间段的baseline作为参考判断event

if iscell(activities)
    rawacti.acq_mean_blocks={};
    rawacti.hab={};
    rawacti.tst={};
    rawacti.acq_1st_block={};
    for ii=1:max(size(activities))
        a=activities{ii};
        for jj=1:size(a,2)
            aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
            rawacti.hab{ii,1}(:,jj)=mean(aa(trial.hab(2):trial.hab(3)-indcut1,:),1);
            rawacti.tst{ii,1}(:,jj)=mean(aa(trial.test(2):trial.test(3)-indcut2,:),1);
            rawacti.acq_1st_block{ii,jj}=aa(trial.acq(2):trial.acq_block_trial:trial.acq(3),:)';
            for zz=1:trial.acq_block_num
                ind1=trial.acq_block_trial*(zz-1)+1+trial.hab(3);
                ind2=trial.acq_block_trial*(zz-1)+1+trial.hab(3)+trial.acq_block_trial-1;
                mean_trial=ind1:ind2;
                if cut_move
                    if ~isempty(intersect(ind_cut_trial,ind1:ind2))
                        inter=intersect(ind_cut_trial,ind1:ind2);
                        mean_trial=setdiff(mean_trial,inter);
                    end
                end
                if ~isempty(mean_trial)
                     rawacti.acq_mean_blocks{ii,jj}(:,zz)=mean(aa(mean_trial,:),1);
                else
                     rawacti.acq_mean_blocks{ii,jj}(:,zz)=[zeros(size(mean(aa(mean_trial,1:frame.us_start-1),1))) mean(aa(:,frame.us_start:end),1)];
                end
            end
        end
    end
    
    area={};
    for ii=1:max(size(activities))
        a=activities{ii};
        for jj=1:size(a,2)
            area.CS_hab_tst{ii,1}(1,jj)=calculate_integrate_dfdf(rawacti.hab{ii,1}(:,jj),area_win_hab_tst,type_cal_area,ref_win);
            area.CS_hab_tst{ii,1}(2,jj)=calculate_integrate_dfdf(rawacti.tst{ii,1}(:,jj),area_win_hab_tst,type_cal_area,ref_win);
            area.CS_acq_block{ii,1}(:,jj)=calculate_integrate_dfdf((rawacti.acq_mean_blocks{ii,jj}),area_win_acq,type_cal_area,ref_win);
            area.CS_acq_1st_block{ii,1}(:,jj)=calculate_integrate_dfdf((rawacti.acq_1st_block{ii,jj}),area_win_acq,type_cal_area,ref_win);
        end
    end
    
elseif ismatrix(activities)
    rawacti.acq_mean_blocks=[];
    rawacti.hab=[];
    rawacti.tst=[];
    rawacti.acq_1st_block=[];
    %按列
    for jj=1:size(activities,2)
        a=activities(:,jj);
            aa=reshape(a,frame.per_cycle,trial.total)';
            rawacti.hab(:,jj)=mean(aa(trial.hab(2):trial.hab(3)-indcut1,:),1);
            rawacti.tst(:,jj)=mean(aa(trial.test(2):trial.test(3)-indcut2,:),1);
            rawacti.acq_1st_block{jj,1}=aa(trial.acq(2):trial.acq_block_trial:trial.acq(3),:)';
            for zz=1:trial.acq_block_num
                ind1=trial.acq_block_trial*(zz-1)+1+trial.hab(3);
                ind2=trial.acq_block_trial*(zz-1)+1+trial.hab(3)+trial.acq_block_trial-1;
                mean_trial=ind1:ind2;
                if cut_move
                    if ~isempty(intersect(ind_cut_trial,ind1:ind2))
                        inter=intersect(ind_cut_trial,ind1:ind2);
                        mean_trial=setdiff(mean_trial,inter);
                    end
                end
                if ~isempty(mean_trial)
                     rawacti.acq_mean_blocks{jj,1}(:,zz)=mean(aa(mean_trial,:),1);
                else
                     rawacti.acq_mean_blocks{jj,1}(:,zz)=[zeros(size(mean(aa(mean_trial,1:frame.us_start-1),1))) mean(aa(:,frame.us_start:end),1)];
                end
            end
    end
        area=[];
    for jj=1:size(activities,2)
            area.CS_hab_tst(:,jj)=calculate_integrate_dfdf([rawacti.hab(:,jj),rawacti.tst(:,jj)],area_win_hab_tst,type_cal_area,ref_win);
            %area.CS_hab_tst(2,jj)=calculate_integrate_dfdf(rawacti.tst(:,jj),area_win_hab_tst,type_cal_area,ref_win);
            area.CS_acq_block(:,jj)=calculate_integrate_dfdf((rawacti.acq_mean_blocks{jj,1}),area_win_acq,type_cal_area,ref_win);
            area.CS_acq_1st_block(:,jj)=calculate_integrate_dfdf((rawacti.acq_1st_block{jj,1}),area_win_acq,type_cal_area,ref_win);
    end
end