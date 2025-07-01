%%
function [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(delta_r_bef1,re_startpoint,fishname,date,outputpath,isbinary)
global fs
global frameb
global trial

% [inputname,inputpath]=uigetfile('G:\data\.mat','behavior_fin');
% load([inputpath,inputname]);
%parameters
frameb.cs_dur = frameb.cs_end-frameb.cs_start;
frameb.csus_interval = frameb.us_start-frameb.cs_start;
%framebnum = trial.total*frameb.per_cycle;
%%
%�ö������ķ������бȽ�
%�ҵ�hab����event(288, 192)
bef_cond_shake={}; bef_cd_shake = [];a=[];b=1;c=1;
for ii = [frameb.cs_dur frameb.csus_interval]
    for jj=1:trial.hab(1)
        a = find(re_startpoint(:,1)==jj);
        for kk = 0:3
            if    find(sum(re_startpoint(a,2)>(frameb.cs_start-ii*(kk-1)-1) & re_startpoint(a,2)<(frameb.cs_end-ii*(kk-1)))) ~= 0;
                if isbinary
                    shake = 1;
                else
                shake = sum(re_startpoint(a,2)>(frameb.cs_start-ii*(kk-1)-1) & re_startpoint(a,2)<(frameb.cs_end-ii*(kk-1)));%1;
                end
            else
                shake=0;
            end
            bef_cd_shake(b,(kk+1)) = shake;
        end
        b=b+1;
        bef_cond_shake{c,1} = bef_cd_shake;
    end
    c=c+1;
end
bef_cond.aftCS_shake = bef_cond_shake{1,1}(:,1);%CS��4.8s
bef_cond.CS_shake = bef_cond_shake{1,1}(:,2);%CS�ڼ�
bef_cond.befCS288_shake = bef_cond_shake{1,1}(:,3);%CSǰ4.8s
bef_cond.befCS576_shake = bef_cond_shake{1,1}(:,4);%CSǰ9.6-4.8s
bef_cond.befCS_shake = (bef_cond.befCS288_shake + bef_cond.befCS576_shake)/2;%ԼΪCSǰȫ��

bef_cond.aftCS192_shake = bef_cond_shake{2,1}((trial.hab(1)+1):(trial.hab(1)*2),1);%ԭΪCSǰ3.2s
bef_cond.CS192_shake = bef_cond_shake{2,1}((trial.hab(1)+1):(trial.hab(1)*2),2);%CS-US�ڼ�
bef_cond.befCS384_shake = bef_cond_shake{2,1}((trial.hab(1)+1):(trial.hab(1)*2),3);
bef_cond.befCS576_shake = bef_cond_shake{2,1}((trial.hab(1)+1):(trial.hab(1)*2),4);
bef_cond.befCS192_shake = (bef_cond.befCS384_shake + bef_cond.befCS576_shake)/2;%ԼΪCSǰȫ��


%�ҵ�acq����event(192)
dur_cd_shake = [];a=[];b=1;
for ii=(trial.hab(1)+1):(trial.hab(1)+trial.acq(1))
    a = find(re_startpoint(:,1)==ii);
    for jj=0:3
        if    find(sum(re_startpoint(a,2)>(frameb.cs_start-frameb.csus_interval*jj-1) & re_startpoint(a,2)<(frameb.us_start-frameb.csus_interval*jj)))~=0;
            if isbinary
                    shake = 1;
            else
                shake = sum(re_startpoint(a,2)>(frameb.cs_start-frameb.csus_interval*jj-1) & re_startpoint(a,2)<(frameb.us_start-frameb.csus_interval*jj));
            end
        else
            shake=0;
        end
        dur_cd_shake(b,(jj+1)) = shake;
    end
    b=b+1;
end

if ~isempty(dur_cd_shake)
    dur_cond.CS192_shake = dur_cd_shake(:,1);%cs-us�ڼ�
    dur_cond.befCS192_shake = dur_cd_shake(:,2);%csǰ4.8S
    dur_cond.befCS384_shake = dur_cd_shake(:,3);%csǰ9.6S
    dur_cond.befCS576_shake = dur_cd_shake(:,4);
    dur_cond.befCS_shake = (dur_cond.befCS192_shake + dur_cond.befCS384_shake)/2;
else
    dur_cond.CS192_shake = [];
    dur_cond.befCS192_shake = [];
    dur_cond.befCS384_shake = [];
    dur_cond.befCS576_shake = [];
    dur_cond.befCS_shake = (dur_cond.befCS192_shake + dur_cond.befCS384_shake)/2;
end

%�ҵ�test����event
aft_cd_shake = [];a=[];b=1;
for ii=(trial.hab(1)+trial.acq(1)+1):trial.total
    a = find(re_startpoint(:,1)==ii);
    for jj=0:3
        if    find(sum(re_startpoint(a,2)>(frameb.cs_start-frameb.cs_dur*(jj-1)-1) & re_startpoint(a,2)<(frameb.cs_end-frameb.cs_dur*(jj-1))))~=0;
             if isbinary
                    shake = 1;
             else
                 shake = sum(re_startpoint(a,2)>(frameb.cs_start-frameb.cs_dur*(jj-1)-1) & re_startpoint(a,2)<(frameb.cs_end-frameb.cs_dur*(jj-1)));
             end
        else
            shake=0;
        end
        aft_cd_shake(b,(jj+1)) = shake;
    end
    b=b+1;
end
if ~isempty( aft_cd_shake)
aft_cond.aftCS_shake = aft_cd_shake(:,1);
aft_cond.CS_shake = aft_cd_shake(:,2);
aft_cond.befCS288_shake = aft_cd_shake(:,3);
aft_cond.befCS576_shake = aft_cd_shake(:,4);
aft_cond.befCS_shake = (aft_cond.befCS288_shake + aft_cond.befCS576_shake)/2;
else
    aft_cond.aftCS_shake = [];
aft_cond.CS_shake = [];
aft_cond.befCS288_shake = [];
aft_cond.befCS576_shake = [];
aft_cond.befCS_shake = (aft_cond.befCS288_shake + aft_cond.befCS576_shake)/2;
end
[p_value dur_cond] = binom_test(re_startpoint,trial,frameb,bef_cond,dur_cond,aft_cond);

%��sepplot����correlation�Ľ�����Ƚ�tail��fin�������?
re_delta_r_bef1=reshape(delta_r_bef1,[frameb.per_cycle,trial.total]);
for ii=1
    h=figure('position',[11,118,1700,860]);
    % set(gcf,'Units','normalized', 'position',get(0,'ScreenSize')) ;
    [h_all, incre_all] = sepplot((1:size(re_delta_r_bef1,1))/fs.behavior,re_delta_r_bef1,repmat([0 0 1],trial.total,1));
    x_range = [0 50 50 0];FontSize=14;
    y_max_all = -ceil(max(-incre_all(:,1)));
    y_incre_mean_all = y_max_all/(size(incre_all,1)-1);
    y_range_hab = [-y_incre_mean_all*0.6 -y_incre_mean_all*0.6 (trial.hab(1)/trial.total)*y_max_all (trial.hab(1)/trial.total)*y_max_all];
    y_range_acq = [(trial.hab(1)/trial.total)*y_max_all (trial.hab(1)/trial.total)*y_max_all ((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all ((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all];
    y_range_test = [((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all ((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all y_max_all+y_incre_mean_all y_max_all+y_incre_mean_all];
    %     imagesc((1:size(re_delta_r_bef1,1))/fs.behavior,1:size(re_delta_r_bef1,2), re_delta_r_bef1', [-0.5 0.5]);
    %     sepplot((1:size(re_delta_r_bef1,1))/fs.behavior, re_delta_r_bef1);
    ylim([y_max_all+y_incre_mean_all -y_incre_mean_all*0.6]);xlim([0 frameb.per_cycle/fs.behavior]);box off
    set(gca,'ytick',[],'ycolor','w','FontSize',20,'LineWidth',3)
    title(['All trials', '    ', date, '-', fishname],'Fontsize',35,'Fontname', 'Times New Roman')
    xlabel('Time(s)','Fontsize',30,'Fontname', 'Times New Roman')
    ylabel('Trial','Fontsize',30,'Fontname', 'Times New Roman')
    hold on
    patch(x_range,y_range_hab,'g','FaceAlpha',.05,'linestyle','none')
    patch(x_range,y_range_acq,'r','FaceAlpha',.1,'linestyle','none')
    patch(x_range,y_range_test,'b','FaceAlpha',.05,'linestyle','none')
    %     plot([1 1],[0 2000],'g','LineStyle','--','LineWidth',1)
    %     plot([1 1],[2000 3000],'b','LineStyle','--','LineWidth',1)
    text(12,(-y_incre_mean_all),'CS','FontSize',20)
    plot([frameb.cs_start/fs.behavior frameb.cs_start/fs.behavior],[0 y_max_all],'k','LineStyle','--','LineWidth',2);
    plot([frameb.cs_end/fs.behavior frameb.cs_end/fs.behavior],[0 y_max_all],'k','LineStyle','--','LineWidth',2);
    plot([frameb.us_start/fs.behavior frameb.us_start/fs.behavior],[(trial.hab(1)/trial.total)*y_max_all ((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all],'r','LineStyle','-','LineWidth',3);
    %     plot([frameb.us_start/fs.behavior frameb.us_start/fs.behavior],[(trial.hab(1)/trial.total)*y_max_all ((trial.hab(1)+10)/trial.total)*y_max_all],'r','LineStyle','-','LineWidth',3);
    %     plot([frameb.us_start/fs.behavior frameb.us_start/fs.behavior],[((trial.hab(1)+12)/trial.total)*y_max_all ((trial.hab(1)+21)/trial.total)*y_max_all],'r','LineStyle','-','LineWidth',3);
    %     plot([frameb.us_start/fs.behavior frameb.us_start/fs.behavior],[((trial.hab(1)+23)/trial.total)*y_max_all ((trial.hab(1)+trial.acq(1))/trial.total)*y_max_all],'r','LineStyle','-','LineWidth',3);
    %�����trial
    yt_all_hab = incre_all(1:2:trial.hab(1),1);
    xt_all_hab = -0.5*ones(size(yt_all_hab));
    txt_all_hab = num2cell(1:2:trial.hab(1));%������Ҫʹ��cell��ʽ���У����ʹ��chr��ʽ�Ļ�������Ҫת�ó�������
    text(xt_all_hab,yt_all_hab,txt_all_hab,'FontSize',FontSize);
    text(-1,y_max_all*trial.hab(1)/(1.02*trial.total),'Habituation','FontSize',25,'rotation',90)
    yt_all_acq = incre_all((trial.hab(1)+1):2:(trial.hab(1)+trial.acq(1)),1);
    xt_all_acq = -0.5*ones(size(yt_all_acq));
    txt_all_acq = num2cell(1:2:trial.acq(1));
    text(xt_all_acq,yt_all_acq,txt_all_acq,'FontSize',FontSize);
    %     yt_all_acq_1 = incre_all((trial.hab(1)+1):6:(trial.hab(1)+trial.acq(1)),1);
    %     xt_all_acq_1 = 0.02*ones(size(yt_all_acq_1));
    %     scatter(xt_all_acq_1,yt_all_acq_1,88,'>','r','filled');
    text(-1,y_max_all*(trial.hab(1)*2+trial.acq(1))/(1.76*trial.total),'Acquisition','FontSize',25,'rotation',90)
    yt_all_test = incre_all((trial.hab(1)+trial.acq(1)+1):2:trial.total,1);
    xt_all_test = -0.5*ones(size(yt_all_test));
    txt_all_test = num2cell(1:2:trial.test(1));
    text(xt_all_test,yt_all_test,txt_all_test,'FontSize',FontSize);
    text(-1,y_max_all*(trial.hab(1)+trial.acq(1)+trial.total)/(1.92*trial.total),'Test','FontSize',25,'rotation',90)
    re_startpoint_sd_incre = re_startpoint;
    index = 0;a=[];
    for ii=1:trial.total
        index = index+1;
        a = sum(re_startpoint_sd_incre(:,1)==ii);
        re_startpoint_sd_incre([index:(index+a-1)],1)=incre_all(ii,1);
        index = index+a-1;
    end
    scatter(re_startpoint_sd_incre(:,2)/fs.behavior,re_startpoint_sd_incre(:,1),[66],'m','filled');
    %����ѧϰǰ���Լ�ѧϰ��spontaneous��CS����������ѧ��
    a=20.2;
    if p_value.CS_15_3 <=0.05 & p_value.aftcond_befCS_CS_3<=0.05
        text(a,(-y_incre_mean_all),'learner(3)','FontSize',FontSize+2)
    else
        text(a,(-y_incre_mean_all),'non-learner(3)','FontSize',FontSize+2)
    end
    if p_value.CS_15_6 <=0.05 & p_value.aftcond_befCS_CS_6<=0.05
        text(a,(-y_incre_mean_all+12*y_incre_mean_all),'learner(6)','FontSize',FontSize+2)
    else
        text(a,(-y_incre_mean_all+12*y_incre_mean_all),'non-learner(6)','FontSize',FontSize+2)
    end
    if 1-p_value.CS_15_3 <=0.05 & 1-p_value.aftcond_befCS_CS_3<=0.05
        text(a,(-y_incre_mean_all+24*y_incre_mean_all),'freezing-learner(3)','FontSize',FontSize+2)
    else
        text(a,(-y_incre_mean_all+24*y_incre_mean_all),'freezing-non-learner(3)','FontSize',FontSize+2)
    end
    if 1-p_value.CS_15_6 <=0.05 & 1-p_value.aftcond_befCS_CS_6<=0.05
        text(a,(-y_incre_mean_all+32*y_incre_mean_all),'freezing-learner(6)','FontSize',FontSize+2)
    else
        text(a,(-y_incre_mean_all+32*y_incre_mean_all),'freezing-non-learner(6)','FontSize',FontSize+2)
    end
    text(a,(-y_incre_mean_all+2*y_incre_mean_all),['p(','cs',') = ',num2str(p_value.CS_15_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+4*y_incre_mean_all),['p(','spon',') = ',num2str(p_value.aftcond_befCS_CS_3)],'FontSize',FontSize)
    %     text(a,(-y_incre_mean_all+6*y_incre_mean_all),['p(','cs-576',') = ',num2str(p_value.aftcond_befCS576_CS_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+6*y_incre_mean_all),['p(','states',') = ',num2str(p_value.states_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+14*y_incre_mean_all),['p(','cs',') = ',num2str(p_value.CS_15_6)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+16*y_incre_mean_all),['p(','spon',') = ',num2str(p_value.aftcond_befCS_CS_6)],'FontSize',FontSize)
    %     text(a,(-y_incre_mean_all+18*y_incre_mean_all),['p(','cs-576',') = ',num2str(p_value.aftcond_befCS576_CS_6)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+18*y_incre_mean_all),['p(','states',') = ',num2str(p_value.states_6)],'FontSize',FontSize)
    %     text(a,(-y_incre_mean_all+26*y_incre_mean_all),['p(','hab-288',') = ',num2str(p_value.befcond_befCS_CS_6')],'FontSize',11)
    text(a,(-y_incre_mean_all+26*y_incre_mean_all),['p(','cs',') = ',num2str(1-p_value.CS_15_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+28*y_incre_mean_all),['p(','spon',') = ',num2str(1-p_value.aftcond_befCS_CS_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+34*y_incre_mean_all),['p(','cs',') = ',num2str(1-p_value.CS_15_6)],'FontSize',FontSize)
    text(a,(-y_incre_mean_all+36*y_incre_mean_all),['p(','spon',') = ',num2str(1-p_value.aftcond_befCS_CS_6)],'FontSize',FontSize)
    %����acq�׶�ÿ��block��ѧϰ֮ǰ�Լ�spontaneous�Ƚϣ���һ��������������(����ļ�����ITI�㹻�ö���ָ�?)
    a=a-5;
    if sum(find(p_value.durcond_befCS192_CS<=0.05))~=0 & sum(find(p_value.dur_learnning_CS_6<=0.05))~=0
        text(a-5,1.14*y_max_all,'preparatory learner(acq-hab(6))','FontSize',10)
    else
        text(a-5,1.14*y_max_all,'non-learner(acq-hab(6))','FontSize',10)
    end
    if sum(find(p_value.durcond_befCS192_CS<=0.05))~=0 & sum(find(p_value.dur_learnning_CS_15<=0.05))~=0
        text(a,1.14*y_max_all,'preparatory learner(acq-hab(all))','FontSize',10)
    else
        text(a,1.14*y_max_all,'non-learner(acq-hab(all))','FontSize',10)
    end
    %����acq�׶Σ�����������block��CS-US�׶ζ��ˣ���û��spontaneousû�У���ô��ѧ��(����ļ�����ITI�������ö���ָ�?)
    re_dur_cond_CS192 = reshape(dur_cond.CS192_shake,trial.acq_block_trial,trial.acq_block_num);
    re_dur_cond_befCS192 = reshape(dur_cond.befCS192_shake,trial.acq_block_trial,trial.acq_block_num);
    shake_csus_loc = find(re_dur_cond_CS192(1,2:trial.acq_block_num)==1);
    if length(shake_csus_loc)>1 & sum(re_dur_cond_befCS192(1,shake_csus_loc))==0
        text(a+5,1.14*y_max_all,'preparatory learner(acq-1)','FontSize',10)
    else
        text(a+5,1.14*y_max_all,'non-learner(acq-1)','FontSize',10)
    end
    text(a,1.08*y_max_all,['p(','acq-192',') = ',num2str(roundn(p_value.durcond_befCS192_CS',-2))],'FontSize',12)
    text(a,1.1*y_max_all,['p(','acq-384',') = ',num2str(roundn(p_value.durcond_befCS384_CS',-2))],'FontSize',12)
    text(a,1.12*y_max_all,['p(','acq-576',') = ',num2str(roundn(p_value.durcond_befCS576_CS',-2))],'FontSize',12)
    hold off
    savefig(h,fullfile(outputpath,[fishname,'behavior statistics-1']));
    %saveas(h,fullfile(outputpath,[fishname,'behavior','.jpeg']));  
    %close(h) 
end
 
%ֻ����hab��test

for ii=1
    h=figure('position',[11,118,1500,600]);
    re_delta_r_bef1_hab_test = re_delta_r_bef1(:,[1:trial.hab(1) (trial.hab(1)+trial.acq(1)+1):trial.total]);
    [h_hab_test, incre_hab_test] = sepplot((1:size(re_delta_r_bef1,1))/fs.behavior, re_delta_r_bef1_hab_test,repmat([0 0 1],trial.total-trial.acq(1),1));
    x_range = [0 50 50 0];
    y_max_hab_test = -ceil(max(-incre_hab_test(:,1)));
    y_incre_mean_hab_test = y_max_hab_test/(size(incre_hab_test,1)-1);
    y_range_hab = [-y_incre_mean_hab_test*0.8 -y_incre_mean_hab_test*0.8 (trial.hab(1)/(trial.hab(1)+trial.test(1)))*y_max_hab_test (trial.hab(1)/(trial.hab(1)+trial.test(1)))*y_max_hab_test];
    y_range_test = [(trial.hab(1)/(trial.hab(1)+trial.test(1)))*y_max_hab_test (trial.hab(1)/(trial.hab(1)+trial.test(1)))*y_max_hab_test y_max_hab_test+y_incre_mean_hab_test y_max_hab_test+y_incre_mean_hab_test];
    %     sepplot((1:size(re_delta_r_bef1,1))/fs.behavior, re_delta_r_bef1_hab_test);
    ylim([y_max_hab_test+y_incre_mean_hab_test -y_incre_mean_hab_test*0.8]);box off
    set(gca,'ytick',[],'ycolor','w','FontSize',20,'LineWidth',3)
    title(['Hab & Test', '    ', date, '-', fishname],'Fontsize',35,'Fontname', 'Times New Roman')
    xlabel('Time(s)','Fontsize',30,'Fontname', 'Times New Roman')
    ylabel('Trial','Fontsize',30,'Fontname', 'Times New Roman')
    hold on
    patch(x_range,y_range_hab,'g','FaceAlpha',.05,'linestyle','none')
    patch(x_range,y_range_test,'b','FaceAlpha',.05,'linestyle','none')
    text(12,(y_max_hab_test*1.02),'CS','FontSize',20)
    plot([frameb.cs_start/fs.behavior frameb.cs_start/fs.behavior],[0 y_max_hab_test],'k','LineStyle','--','LineWidth',2);
    plot([frameb.cs_end/fs.behavior frameb.cs_end/fs.behavior],[0 y_max_hab_test],'k','LineStyle','--','LineWidth',2);
    %�����trial
    yt_hab = incre_hab_test(1:2:trial.hab(1),1);
    xt_hab = -0.5*ones(size(yt_hab));
    txt_hab = num2cell(1:2:trial.hab(1));%������Ҫʹ��cell��ʽ���У����ʹ��chr��ʽ�Ļ�������Ҫת�ó�������
    text(xt_hab,yt_hab,txt_hab,'FontSize',FontSize)
    text(-1,y_max_hab_test*trial.hab(1)/(1.6*(trial.hab(1)+trial.test(1))),'Habituation','FontSize',25,'rotation',90)
    yt_test = incre_hab_test((trial.hab(1)+1):2:end,1);
    xt_test = -0.5*ones(size(yt_test));
    txt_test = num2cell(1:2:trial.test(1));%������Ҫʹ��cell��ʽ���У����ʹ��chr��ʽ�Ļ�������Ҫת�ó�������
    text(xt_test,yt_test,txt_test,'FontSize',FontSize)
    text(-1,y_max_hab_test*(trial.hab(1)*2+trial.test(1))/(1.9*(trial.hab(1)+trial.test(1))),'Test','FontSize',25,'rotation',90)
    loc_re_startpoint_sd_hab_test = find(re_startpoint(:,1)<(trial.hab(1)+1)|re_startpoint(:,1)>(trial.hab(1)+trial.acq(1)));
    re_startpoint_sd_incre_hab_test = re_startpoint(loc_re_startpoint_sd_hab_test,:);
    index = 0;
    for ii=1:trial.hab(1)
        index = index+1;
        a = sum(re_startpoint_sd_incre_hab_test(:,1)==ii);
        re_startpoint_sd_incre_hab_test([index:(index+a-1)],1)=incre_hab_test(ii,1);
        index = index+a-1;
    end
    for ii=(trial.hab(1)+trial.acq(1)+1):trial.total
        index = index+1;
        a = sum(re_startpoint_sd_incre_hab_test(:,1)==ii);
        re_startpoint_sd_incre_hab_test([index:(index+a-1)],1)=incre_hab_test((ii-trial.acq(1)),1);
        index = index+a-1;
    end
    scatter(re_startpoint_sd_incre_hab_test(:,2)/fs.behavior,re_startpoint_sd_incre_hab_test(:,1),[66],'m','filled');
    a=30.2;
    if p_value.CS_15_3 <=0.05 & p_value.aftcond_befCS_CS_3<=0.05
        text(a,(-y_incre_mean_hab_test),'learner(3)','FontSize',20)
    else
        text(a,(-y_incre_mean_hab_test),'non-learner(3)','FontSize',20)
    end
    if p_value.CS_15_6 <=0.05 & p_value.aftcond_befCS_CS_6<=0.05
        text(a,(-y_incre_mean_hab_test+6*y_incre_mean_hab_test),'learner(6)','FontSize',20)
    else
        text(a,(-y_incre_mean_hab_test+6*y_incre_mean_hab_test),'non-learner(6)','FontSize',20)
    end
    xlim([0 frameb.per_cycle/fs.behavior])
    text(a,(-y_incre_mean_hab_test+1*y_incre_mean_hab_test),['p(','cs',') = ',num2str(p_value.CS_15_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_hab_test+2*y_incre_mean_hab_test),['p(','spon',') = ',num2str(p_value.aftcond_befCS_CS_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_hab_test+3*y_incre_mean_hab_test),['p(','states',') = ',num2str(p_value.states_3)],'FontSize',FontSize)
    text(a,(-y_incre_mean_hab_test+7*y_incre_mean_hab_test),['p(','cs',') = ',num2str(p_value.CS_15_6)],'FontSize',FontSize)
    text(a,(-y_incre_mean_hab_test+8*y_incre_mean_hab_test),['p(','spon',') = ',num2str(p_value.aftcond_befCS_CS_6)],'FontSize',FontSize)
    text(a,(-y_incre_mean_hab_test+9*y_incre_mean_hab_test),['p(','states',') = ',num2str(p_value.states_6)],'FontSize',FontSize)
    hold off
    %savefig(h,fullfile(outputpath,[fishname,'behavior statistics-2']));
    close(h)
end
