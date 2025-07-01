%add non-learner
%processing results
%add n
%前后比较，learner/non-learner比较
%
%% load
clc;clear all;
actbatch={};envbatch={};behavbatch={};parabatch={};
actbatch=batch_code('txt of <activties_new> path') ;
envbatch=batch_code('txt of <location> path') ;
behavbatch=batch_code('txt of <behavior> path') ;
parabatch=batch_code('txt of <para> path') ;
savebatch=batch_code('txt of <save> path') ;
savepath_all=uigetdir('G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\','savepath');
%save([savepath_all '\path'],'actbatch','envbatch','behavbatch','parabatch','-v7.3');
global collength
collength=2100;
%correct actpath
name='activities_aft_process';
for tt=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{tt});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{tt}=fullfile(fpath,[fname,ext]);
    end
end
K=[2,3,4,6,9,12,13];
[act_all,act_all_hab,act_all_acq,act_all_acq_block_mean,act_all_tst,env_all,index_all]=Get_all_act(K,collength,true,true,actbatch,envbatch,behavbatch,parabatch,savepath_all);

for ss=1
    ind=find(env_all.supervoxel(:,3)==1);
    act_all(ind,:)=[];
    act_all_hab(:,:,ind)=[];act_all_tst(:,:,ind)=[];act_all_acq(:,:,:,ind)=[];act_all_acq_block_mean(:,:,ind)=[];
    env_all.supervoxel(ind,:)=[];
    index_all(ind)=[];
end
for ss=[9]
    ind=find(index_all==ss);
    act_all(ind,:)=[];
    act_all_hab(:,:,ind)=[];act_all_tst(:,:,ind)=[];act_all_acq(:,:,:,ind)=[];act_all_acq_block_mean(:,:,ind)=[];
    env_all.supervoxel(ind,:)=[];
    index_all(ind)=[];
end
act_all_hab_mean=mean(act_all_hab,2);act_all_tst_mean=mean(act_all_tst,2);
unique(index_all)
%% plot para
colorCS=[0.5 0.5 0.9];
%% regressor
%[CS*7/US*5/motor*2/autoclus]
%CS: HAB/TRIANING/TEST
%MOV:SPON/INDUCED/
%STATE:
%US
load(parabatch{1});num_fish_freq=2;
stimcorr.CS=zeros(size(act_all,1),trial.acq_block_num+2);
stAvrcorr.CS=zeros(size(act_all,1),trial.acq_block_num+2);
stimcorr.CSCS=zeros(size(act_all,1),trial.acq_block_num+2);
stAvrcorr.CSCS=zeros(size(act_all,1),trial.acq_block_num+2);

stimcorr.US=zeros(size(act_all,1),trial.acq_block_num);
motorcorr.all=zeros(size(act_all,1),1);
motorcorr.CS=zeros(size(act_all,1),1);
motorcorr.spon=zeros(size(act_all,1),1);
motorcorr.sep_all=zeros(size(act_all,1),trial.acq_block_num+2);
motorcorr.sep_CS=zeros(size(act_all,1),trial.acq_block_num+2);
motorcorr.sep_spon=zeros(size(act_all,1),trial.acq_block_num+2);

reg_para.prct_const=2;reg_para.reg_thres=0.1;reg_para.reg_thres2=0.6;%para
for kk=unique(index_all)'
    beh=load(behavbatch{kk});load(parabatch{kk});
    if ~isfield(beh,'re_startpoint_sd')
        beh.re_startpoint_sd=beh.re_start_end_point;
    end
    ind_in_this_fish=find(index_all==kk);
    for ss=1 % bulid stim regressor
        stimCS=zeros(1,frame.per_cycle*trial.acq_block_trial);
        for tt=1:size(stimCS,2)/frame.per_cycle
            stimCS((tt-1)*frame.per_cycle+frame.cs_start:(tt-1)*frame.per_cycle+frame.cs_end)=1;%23
        end
        stimUS=zeros(1,frame.per_cycle*trial.acq_block_trial);
        for tt=1:size(stimUS,2)/frame.per_cycle
            stimUS((tt-1)*frame.per_cycle+frame.us_start)=3;%23
        end
        stim=zeros(5,trial.acq_block_trial*frame.per_cycle);stim(1,:)=stimCS;stim(2,:)=stimUS;
        [stregressors,~,~,~] = GetMotorRegressor(stim,num_fish_freq);regress.CS_acq_block=stregressors(1,1).im;regress.US_acq_block=stregressors(1,2).im;
        figure,plot(stimCS);hold on; plot(regress.CS_acq_block,'r');hold on;plot(stimUS);hold on; plot(regress.US_acq_block,'r');
        %corr(regress.CS_acq_block,regress.US_acq_block')
    end
    stimcorr_k=nan(length(ind_in_this_fish),7);stAvrcorr_k=nan(length(ind_in_this_fish),7);stimcorr_ku=nan(length(ind_in_this_fish),7);
    stimcorr_k_cs=nan(length(ind_in_this_fish),7);stAvrcorr_k_cs=nan(length(ind_in_this_fish),7);
    motorcorr_k_all=nan(length(ind_in_this_fish),7);motorcorr_k_CS=nan(length(ind_in_this_fish),7);motorcorr_k_spon=nan(length(ind_in_this_fish),7);
    behavior=getmulti_motor_regresspr(beh.re_startpoint_sd,beh.startpoint_sd,beh.y_3sd,frameb,fs);
    [mregressors,~,~,~] = GetMotorRegressor(behavior,6);
    regress.motor_all=mregressors(1,1).im;regress.motor_CS=mregressors(1,2).im;regress.motor_spon=mregressors(1,3).im;
    time_be=0:1/fs.behavior:(length(beh.y_3sd)-1)/fs.behavior;
    time_ca=0:fs.ca:(size(act_all,2)-1)*fs.ca;
    regress.motor_all_down = interp1(time_be,regress.motor_all,time_ca);regress.motor_all_down=regress.motor_all_down(1:collength);
    regress.motor_CS_down = interp1(time_be,regress.motor_CS,time_ca);regress.motor_CS_down=regress.motor_CS_down(1:collength);
    regress.motor_spon_down = interp1(time_be,regress.motor_spon,time_ca);regress.motor_spon_down=regress.motor_spon_down(1:collength);
    figure,plot(time_ca,regress.motor_all_down ,'k');hold on;plot(time_ca,regress.motor_CS_down ,'r');hold on;plot(time_ca,regress.motor_spon_down ,'b');hold on;
    M=act_all(ind_in_this_fish,:);
    figure, plot_rawtrace(M',[1,collength],[],fs,frame,trial,startpoint);title(['Fish ',num2str(kk)]);
    for ss=1 %motor
    [motorcorr.all(ind_in_this_fish,:),~,~]=get_regressor_M_corr(M, regress.motor_all_down,reg_para,frame,1:frame.per_cycle);
    [motorcorr.CS(ind_in_this_fish,:),~,~]=get_regressor_M_corr(M, regress.motor_CS_down,reg_para,frame,1:frame.per_cycle);
    [motorcorr.spon(ind_in_this_fish,:),~,~]=get_regressor_M_corr(M,regress.motor_spon_down,reg_para,frame,1:frame.per_cycle);
    end
    for tt=1:trial.acq_block_num+2
        if tt>=2 && tt<=6
            frame_ind_cs=1:frame.us_start-1;isus=true;
        else
            frame_ind_cs=1:frame.per_cycle;isus=false;
        end
        switch tt
            case 1
                trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
            case {2,3,4,5,6}
                trial_ind=trial.acq(2)+trial.acq_block_trial*(tt-2):trial.acq(2)+trial.acq_block_trial*(tt-1)-1;
            case 7
                trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
        end
        M=act_all(ind_in_this_fish,(trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end)*frame.per_cycle));
        M_0=M;
        stim_output=[regress.CS_acq_block];%stim_output=reshape(stim_output',frame.per_cycle,[]);stim_output=stim_output(frame_ind_cs,:);stim_output=reshape(stim_output,[],1)';
        [stimcorr_k(:,tt),stAvrcorr_k(:,tt),~]=get_regressor_M_corr( M_0,stim_output,reg_para,frame,size(M_0,2)/length(trial_ind));
        M_0=reshape(M',frame.per_cycle,length(trial_ind),[]); M_0= M_0(frame_ind_cs,:,:); M_0=reshape( M_0,length(frame_ind_cs)*length(trial_ind),[])';
        stim_output=[regress.CS_acq_block];stim_output=reshape(stim_output',frame.per_cycle,[]);stim_output=stim_output(frame_ind_cs,:);stim_output=reshape(stim_output,[],1)';
        [stimcorr_k_cs(:,tt),stAvrcorr_k_cs(:,tt),~]=get_regressor_M_corr( M_0,stim_output,reg_para,frame,size(M_0,2)/length(trial_ind));

        [stimcorr_ku(:,tt),~,~]=get_regressor_M_corr(M,regress.US_acq_block,reg_para,frame,1:frame.per_cycle);
        motor_output=regress.motor_all_down((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end)*frame.per_cycle));
        [motorcorr_k_all(:,tt),~,~]=get_regressor_M_corr(M, motor_output,reg_para,frame,1:frame.per_cycle);
        motor_output=regress.motor_CS_down((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end)*frame.per_cycle));
        [motorcorr_k_CS(:,tt),~,~]=get_regressor_M_corr(M, motor_output ,reg_para,frame,1:frame.per_cycle);
        motor_output=regress.motor_spon_down((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end)*frame.per_cycle));
        [motorcorr_k_spon(:,tt),~,~]=get_regressor_M_corr(M,motor_output,reg_para,frame,1:frame.per_cycle);
    end %%ii
    %figure,hist(stAvrcorr_k(:,1));hold on;hist(stAvrcorr_k(:,7));
    stimcorr.CS(ind_in_this_fish,:)=stimcorr_k;
    stAvrcorr.CS(ind_in_this_fish,:)=stAvrcorr_k;
    stimcorr.CSCS(ind_in_this_fish,:)=stimcorr_k_cs;
    stAvrcorr.CSCS(ind_in_this_fish,:)=stAvrcorr_k_cs;
    stimcorr.US(ind_in_this_fish,:)=stimcorr_ku(:,2:end-1);
    motorcorr.sep_all(ind_in_this_fish,:)=motorcorr_k_all;
    motorcorr.sep_CS(ind_in_this_fish,:)=motorcorr_k_CS;
    motorcorr.sep_spon(ind_in_this_fish,:)=motorcorr_k_spon;
end

%regressor_profile=cat(2,stimcorr.CS,stimcorr.US,motorcorr.all,motorcorr.CS,motorcorr.spon,motorcorr.sep_all,motorcorr.sep_CS,motorcorr.sep_spon);%all
%regressor_profile=cat(2,stimcorr.CS,stimcorr.US,stAvrcorr.CSCS);%only stim
regressor_profile=cat(2,stimcorr.CS(:,[1,7]),stAvrcorr.CS(:,[1,7]));%only stim
regressor_profile=cat(2,stimcorr.CSCS,stAvrcorr.CSCS);%only stim
%regressor_profile=stAvrcorr.CSCS;

regressor_profile_n=normalize(regressor_profile,1,'zscore');figure,imagesc(regressor_profile_n ,[-1 1]);colorbar;

[i,j]=ind2sub(size(regressor_profile_n),find(isnan(regressor_profile_n)));
regressor_profile_cut_nan=regressor_profile_n;regressor_profile_cut_nan(:,unique(j))=[];
[coeff,score,~,~,explained,mu] =  pca(regressor_profile_cut_nan,'Rows','complete');
sum_explained = 0;
idx = 0;
while sum_explained < 95
    idx = idx + 1;
    sum_explained = sum_explained + explained(idx);
end
idx
regressor_profile_n_p = score(:,1:idx);
figure,imagesc(regressor_profile_n_p ,[-1 1]);colorbar;
save([savepath_all '\regressor_profile.mat'],'stimcorr','stAvrcorr','motorcorr','reg_para','regressor_profile','regressor_profile_n_p');
%% k_means
for ii=40
    ii
[idx_kmeans,~,sumd,~] = kmeans(regressor_profile_n_p,ii,'Distance','correlation','Replicates',10,'Display','final','MaxIter',1000);
end
cIX=[1:length(idx_kmeans)];
save([savepath_all,'\idx_kmeans_',num2str(length(unique(idx_kmeans)),'%02d'),'.mat'],'idx_kmeans');
%idx_kmeans(find(regressor_profile(:,1)>=0.2))=1;idx_kmeans(find(regressor_profile(:,1)<0.2))=2;
%[idx_kmeans,~,numK,~]=find_best_k_in_range(regressor_profile_n_p,19:20);
cIX=autocluster_raw.cIX;
idx_kmeans=autocluster_raw.gIX;
%% sort gIX
issort=true;ismerge=false;
if issort
    A=act_all_acq_block_mean;cat(2,act_all_hab_mean,act_all_tst_mean);
    m=[];k=[];
    for kk=unique(idx_kmeans)'
        ind=cIX(find(idx_kmeans==kk));
        a=mean(A(:,:,ind),3);
        m(:,kk)=mean(a(frame.cs_start:frame.us_start,:));
        p = polyfit([1:size(A,2)]',m(:,kk),1);
        k(:,kk)=p(1);
    end
    [K,I]=sort(k);
    gIX_s=[];cIX_s=[];kk=1;
    for ii=I
        ind=cIX(find(idx_kmeans==ii));
        gIX_s=[gIX_s;kk*ones(size(ind))];
        cIX_s=[cIX_s;ind];
        kk=kk+1;
    end
else
    gIX_s=idx_kmeans;
    cIX_s=cIX;
end
% cIX(find(autocluster_raw.gIX==2));gIX_s((find(cIX_s==ans(300))))
% idx_kmeans(1000)
% gIX_s(find(cIX_s==1000))
%% merge gIX
if ismerge
   load('G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_nonlearner_cutmove_20190921\autoregress_k41_onlystimCSCS_pca_merge_matrix.mat');%行：每个merge_ind;列:merge后个数
    gIX_m=nan(size(gIX_s));
    for ii=1:size(merge_matrix,1)
        merge_ind=merge_matrix(ii,:);
        for jj=1:length(merge_ind)
            ind=find(gIX_s==merge_ind(jj));
            if ~isempty(ind)
                gIX_m(ind)=ii;
            end
        end
    end
else
    merge_matrix=[];gIX_m=gIX_s;
end
%% correct
for ii=[23,27,18,28,29,13,14,21,11,26,4,6,16,25,35,8,10,12]
cIX_s(find(gIX_m==ii))=[];
gIX_m(find(gIX_m==ii))=[];
end
[K,I]=sort(unique(gIX_m)');
gIX_c=[];cIX_c=[];kk=1;
for ii=K
    ind=cIX_s(find(gIX_m==ii));
    gIX_c=[gIX_c;kk*ones(size(ind))];
    cIX_c=[cIX_c;ind];
    kk=kk+1;
end
cIX_s=cIX_c;gIX_m=gIX_c;
unique(gIX_m)
%% plot para
gIX=gIX_m;
cIX=cIX_s;
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX));
ind_control=setdiff(1:length(index_all),cIX)';%length(ind)+length(cIX)== length(index_all) 
clr_control=[0,0,0];
cIX_all=[cIX;ind_control];
gIX_all=[gIX;(max(gIX)+1)*ones(size(ind_control))];
clrmap_all=cat(1,clrmap,clr_control);
length(unique(gIX))
%% build seed
for ii=1
[~,~,trace_for_ttest_hab]=plot_test_2(act_all_hab_mean,cIX,gIX,frame,clrmap_iii,2,[-0.02 0.04],false,true);
[~,~,trace_for_ttest_tst]=plot_test_2(act_all_tst_mean,cIX,gIX,frame,clrmap_iii,2,[-0.02 0.04],false,true);
[~,~,trace_for_ttest_acq]=plot_test_2(act_all_acq_block_mean,cIX,gIX,frame,clrmap_iii,3,[-0.005 0.01],true,true);
trace_seed=[];
for ii=unique(gIX)'
    a=mean(trace_for_ttest_hab{ii},1)';
    b=mean(trace_for_ttest_tst{ii},1)';
    c=mean(trace_for_ttest_acq{ii},3);
    trace_seed(:,:,ii)=cat(2,a,c,b);
    %figure,plot(trace_seed(:,:,2));
end
corr_with_seed=nan(size(A,3),length(unique(gIX)));
for ii=unique(gIX)'
    a=cat(2,act_all_hab_mean,act_all_acq_block_mean,act_all_tst_mean);a=reshape(a,frame.per_cycle*7,[]);
    b=reshape(trace_seed(:,:,ii),[],1);
    corr_with_seed(:,ii)=corr(a,b);
end
corr_with_seed_n=normalize(corr_with_seed,1);
[idx_kmeans_seed,~,sumD,~] = kmeans(corr_with_seed_n,50,'Distance','correlation','Replicates',1,'Display','final','MaxIter',1000);
%% autocluster
addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
M_0=reshape(act_all_acq_block_mean(frame_ind,:,:),length(frame_ind)*size(act_all_acq_block_mean,2),[])';
M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
M_0=regressor_profile_n_p;
%%%clustering
numK=20;
%numK=max(idx_kmeans_seed);
% ind=index_all;%find(index_all==2);
% ind=ind([1:20:length(ind)]);
% [gIXx,~,numK,~]=find_best_k_in_range(M_0(ind,:),3:15);
% length(unique(gIXx))
% cIXx=ind;cIX_reg = (1:size(M_0,1))';%ind;
if isempty(numK)
    numK=20;
end
masterthres=0.5;
para=struct;
para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',80);
%cIX=(1:size(M_0,1))';gIX=(1:size(M_0,1))'; cIX_reg = (1:size(M_0,1))';
cIXx=1:size(M_0,1);gIXx=idx_kmeans;cIX_reg = (1:size(M_0,1))';
for tt=1:5
%[cIX,gIX] = AutoClustering(randi(size(M_0,1),[1,50000]),[1:50000],M_0,cIX_reg,1,para,1,0.7);
[cIX,gIX] = AutoClustering(cIXx,gIXx,M_0,cIX_reg,0,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
length(unique(gIX))
ind=[1:length(cIX)];%randi(length(cIX),[floor(length(cIX)/5),1]);
cIXx=cIX(ind);
gIXx=gIX(ind);
end
length(unique(gIX))
autocluster_raw.cIX=cIX;autocluster_raw.gIX=gIX;%!!!!!!!!!!!!!!!!!!!!!!!!
autocluster_raw.index_all=index_all;%!!!!勿覆盖
end
%% plot
%cIX(35670:35679)=[];gIX(35670:35679)=[];
vi='on';'off';
for iii=1
    switch iii
        case 1
            cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
        case 2
            cIX_iii=cIX_all;gIX_iii=gIX_all;clrmap_iii=clrmap_all;
    end
    outputpath=checkpath([savepath_all '\figures_regressor_profile' '_' num2str(iii,'%02d')]);
    %count
    figure,h=histogram(gIX_iii,'binedge',[1:max(unique(gIX_iii))+1]);h=bar(h.Values/length(gIX_all));h.FaceColor = 'flat';
    for ii=unique(gIX_iii)
        h.CData(ii,:) = clrmap_iii(ii,:);
    end
    xlim([0.5 max(unique(gIX_iii))+0.5]);xlabel('# Cluster');ylabel('Fraction');set(gca,'fontsize',15,'xcolor','w','ycolor','w')
    %1
    [h,gIX_2, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,regressor_profile_n_p(cIX_iii,:));saveas(h,[outputpath '\f1'],'fig');
    figure,imagesc(regressor_profile_n_p(cIX_iii,:) ,[-2 2]);colorbar;colormap('parula');
    %2
    a=gIX_iii;a(find(gIX_iii==13))=[];b=cIX_iii;b(find(gIX_iii==13))=[];
     B=act_all_acq(frame.cs_start:frame.us_start-1,:,:,:);B=reshape(B,length([frame.cs_start:frame.us_start-1])*trial.acq(1),[])';
     [h]=pushbutton_popupplot_Callback(B,b,a,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[outputpath '\f2'],'fig');
     B=act_all_acq_block_mean(frame.cs_start:frame.us_start-2,:,:,:);B=reshape(B,length([frame.cs_start:frame.us_start-2])*trial.acq_block_num,[])';
     [h]=pushbutton_popupplot_Callback(B,b,a,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[outputpath '\f2'],'fig');
    %3
%     h=figure;B=cat(1,area_all.CS_hab_tst(1,:),area_all.CS_acq_block,area_all.CS_hab_tst(2,:));%B=[area.CS_hab_tst(1,:); area.CS_acq_block ;area.CS_hab_tst(2,:)];
%     a=[-0.01 1];plot_test_1([1:size(B,1)],B,cIX_iii,gIX_iii,clrmap_iii,a,frame,fs.ca,2);saveas(h,[outputpath '\f3'],'fig');
    %4
    %[h,ratio,~]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,1,[-0.02 0.04],true,true);saveas(h,[outputpath '\f41'],'fig');
    [h,ratio,~]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.02 0.04],true,true);saveas(h,[outputpath '\f42'],'fig');
    [h,ratio,trace_for_ttest_acq]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,3,[-0.005 0.02],true,true);saveas(h,[outputpath '\f42'],'fig');
    %5
    [h,~,trace_for_ttest_hab]=plot_test_2(act_all_hab_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.01 0.02],false,true);saveas(h,[outputpath '\f5'],'fig');
    %6
    [h,~,trace_for_ttest_test]=plot_test_2(act_all_tst_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.01 0.02],false,true);saveas(h,[outputpath '\f6'],'fig');
    %
    plot_test_4(act_all_hab_mean,act_all_tst_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.01 0.02],false,true);
    %7
    [h1,h2,ratio_fish]=plot_test_3(act_all',cIX_iii,gIX_iii,clrmap_iii,index_all,envbatch,env_all,outputpath,fs,stimCS,stimUS);
end
%% ttest
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
p=nan(1,length(unique(gIX_iii)));h_test_clus=nan(1,length(unique(gIX_iii)))';
for kk=unique(gIX_iii)'
    x=mean(trace_for_ttest_hab{kk}(:,frame.cs_start:frame.cs_end),2);
    y=mean(trace_for_ttest_test{kk}(:,frame.cs_start:frame.cs_end),2);
    [p(kk),h_test_clus(kk) ]= ranksum(x,y,'tail','left');
    %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
end
p=nan(1,length(unique(gIX_iii)));h_hab_base_clus=nan(1,length(unique(gIX_iii)))';
for kk=unique(gIX_iii)'
    y=mean(trace_for_ttest_hab{kk}(:,frame.cs_start:frame.cs_end),2);
    x=mean(trace_for_ttest_hab{kk}(:,1:frame.cs_start-1),2);
    [p(kk),h_hab_base_clus(kk) ]= ranksum(x,y,'tail','left');
    %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
end
h_acq_clus=nan(length(unique(gIX_iii)),size(trace_for_ttest_acq{1},2));
for kk=unique(gIX_iii)'
    x=mean(trace_for_ttest_hab{kk}(:,frame.cs_start:frame.us_start-2),2);
    for ii=1:size(trace_for_ttest_acq{1},2)
    y=squeeze(mean(trace_for_ttest_acq{kk}(frame.cs_start:frame.us_start-2,ii,:),1));
    [~,h_acq_clus(kk,ii) ]= ranksum(x,y,'tail','left');
    end
    %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
end
%% count fraction per region
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap_all;
num_clust=unique(gIX_iii);ii=1;
fraction_in_region_in_clust=[];num_in_region_in_clust_all=[];
for kk=unique(index_all)'
    ind_in_this_fish=(find(index_all==kk));
    ind_in_this_clus=1:length(gIX_iii);%
    cIX_in_this_clus=cIX_iii(ind_in_this_clus);%
    [IX,i_ind_in_this_fish,i_cIX]=intersect(ind_in_this_fish' ,cIX_in_this_clus');
    %unique(index_all(IX))
    unique(index_all(cIX_iii(ind_in_this_clus(i_cIX))));
    unique(gIX_iii(ind_in_this_clus(i_cIX)));
    %e=load(envbatch{kk});
    if exist([savebatch{kk} '\brain arearegion_mask.mat'])
    load([savebatch{kk} '\brain arearegion_mask.mat']);[savebatch{kk} '\brain arearegion_mask.mat']
    %figure,scatter( env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2), env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1))
%     [~,fraction_in_region_in_clust(:,:,ii),~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),[1:length(i_cIX)],...
%         [env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),3)],clrmap_iii(num_clust,:),true);
    [~,fraction_in_region_in_clust(:,:,ii),~,~,num_in_region_in_clust]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),cIX_iii(ind_in_this_clus(i_cIX)),...
       [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
       [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap_iii(num_clust,:),true);
   num_in_region_in_clust_all(:,:,ii)=num_in_region_in_clust{1};
   title(['Fish ' num2str(kk,'%02d')]);
    ii=ii+1;
    end
end
%% difference of fraction between learner and non-learner
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
l_fraction_in_region_in_clust=mean(fraction_in_region_in_clust(:,:,[1,2,3,4]),3,'omitnan');sum(l_fraction_in_region_in_clust,2,'omitnan')
nl_fraction_in_region_in_clust=mean(fraction_in_region_in_clust(:,:,5:end),3,'omitnan');
%d_fraction_in_region_in_clust=l_fraction_in_region_in_clust;
%d_fraction_in_region_in_clust=[l_fraction_in_region_in_clust-nl_fraction_in_region_in_clust];
d_fraction_in_region_in_clust=[l_fraction_in_region_in_clust-nl_fraction_in_region_in_clust]./(l_fraction_in_region_in_clust+nl_fraction_in_region_in_clust);
a=d_fraction_in_region_in_clust(:,[4,10,20,21,23]);%a(find(a<0.5))=0;
%a=l_fraction_in_region_in_clust(:,[4,10,21,23,20]);%a(find(a<0.5))=0;
figure,imagesc(a');colorbar;colormap('jet');set(gca,'clim',[-1 1],'xtick',[1:21],'ylim',[0.5 length(a)+0.5],'xticklabel',Lable,'XTickLabelRotation',45)
min(d_fraction_in_region_in_clust(:))
max(d_fraction_in_region_in_clust(:))
a=squeeze(sum(fraction_in_region_in_clust,1,'omitnan'));figure,h=bar([1:17],a(1:17,:),'stacked','xdata',[1:17]);
%count
Lable={'l pallium rostral'; 'r pallium rostral';'l pallium middle';'r pallium middle';'l pallium lateral';'r pallium lateral';...
    'l habenula';'r habenula';...
    'l tectum rostral';'r tectum rostral';'l tectum middle'; 'r tectum middle';'l tectum caudal';'r tectum caudal';...
    'l tegmentum';'r tegmentum';...
    'l cerebellum';'r cerebellum';...
    'l hindbrain';'r hindbrain';...
    'longitude'
    };c = categorical(Lable)';
figure('name','difference of fraction between learner and non-learner'),
for ii=1:size(d_fraction_in_region_in_clust,2)
    subplot(ceil(size(d_fraction_in_region_in_clust,2)/4),4,ii);
    if ii<size(d_fraction_in_region_in_clust,2)
        h=bar(d_fraction_in_region_in_clust(:,ii),'FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);set(gca,'xticklabel',[]);    
         set(gca,'xtick',[1:1:length(reg_name)]);
         xlim([0.5 length(reg_name)+0.5]);grid on;
    else
        h=bar(c,d_fraction_in_region_in_clust(:,ii),'FaceColor','flat','XDataMode','manual','xdata',[1:length(reg_name)]);
       set(gca,'xticklabels',Lable,'XTickLabelRotation',45);grid on;
    end
    h.CData = clrmap_iii(ii,:);
    title([num2str(ii,'%02d')]);
    ylim([-1 1]);
end
%% different between learner and nonlearner
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
num_l=find(index_all==1 | index_all==2 |index_all==4);
num_nl=find(index_all==13);
bl_fraction=nan(1,length(unique(gIX_iii)));bn_fraction=nan(1,length(unique(gIX_iii)));
for ii=unique(gIX_iii)'
    a=cIX_iii(find(gIX_iii==ii));
    b=intersect(a,num_l); bl_fraction(ii)=length(b)/length(num_l);
    b=intersect(a,num_nl); bn_fraction(ii)=length(b)/length(num_nl);
end
figure,h=bar([bl_fraction;bn_fraction]','FaceColor','flat');set(gca,'xtick',[1:1:max(unique(gIX_iii))]);
h(1).CData=[1 0 0];h(2).CData=[0 0 1];
xlabel('# Cluster');ylabel('Fraction');
figure,h=bar((bl_fraction-bn_fraction)./(bl_fraction+bn_fraction),'FaceColor','flat');set(gca,'xtick',[1:1:max(unique(gIX_iii))]);
for ii=unique(gIX_iii)'
    h.CData(ii,:)=clrmap_iii(ii,:);
end
xlabel('# Cluster');ylabel('Diff. Fraction');ylim([-1 1]);
%% count special cluster fraction per region
cIX_iii=cIX;gIX_iii=gIX;
clrmap_iii=repmat(clrmap(20,:),length(unique(gIX_iii)),1);num_clust=[17];
for kk=unique(index_all)'
    ind_in_this_fish=(find(index_all==kk));
    ind_in_this_clus=find(gIX_iii== 13);%1:length(gIX_iii);%
    cIX_in_this_clus=cIX_iii(ind_in_this_clus);%
    [IX,i_ind_in_this_fish,i_cIX]=intersect(ind_in_this_fish' ,cIX_in_this_clus');
    %unique(index_all(IX))
%     unique(index_all(cIX_iii(ind_in_this_clus(i_cIX))))
%     unique(gIX_iii(ind_in_this_clus(i_cIX)))
    %e=load(envbatch{kk});
    load([savebatch{kk} '\brain arearegion_mask.mat']);
    [~,~,~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),cIX_iii(ind_in_this_clus(i_cIX)),...
        [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
        [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap_iii(num_clust,:),true);
    title(['Fish ' num2str(kk,'%02d')]);
end
%% motor trigger avg
clrmap_iii=clrmap;cIX_iii=cIX;gIX_iii=gIX;num_clust=unique(gIX_iii)';
color={};color{12}=[0 0 0];color{13}=[0.5 0.5 0.5];
color{2}=[1 0 0];color{3}=[0 0 1];color{4}=[1 0 1];color{6}=[0.2 1 0.2];
for ss=1:2
    a=[];
    figure,
    zz=1;
    for ii=num_clust
        subplot(ceil(length(num_clust)/4),4,zz), line([6 6],[-0.02 0.03]);hold on;
        for kk=unique(index_all)'
            ind_env_in_this_fish=find(index_all==kk);ind_cls=cIX_iii(find(gIX_iii==ii));
            [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(ind_env_in_this_fish ,ind_cls);
            beh=load(behavbatch{kk});
            if ~isfield(beh,'re_startpoint_sd')
                beh.re_startpoint_sd=beh.re_start_end_point;
            end
            switch ss %before
                case 1
                    ind=find(beh.re_startpoint_sd(:,1)<=trial.hab(3) & beh.re_startpoint_sd(:,1)>=trial.hab(2) & beh.re_startpoint_sd(:,2)<frameb.cs_start);
                case 2 %after
                    ind=find(beh.re_startpoint_sd(:,1)<=trial.test(3) & beh.re_startpoint_sd(:,1)>=trial.test(2)& beh.re_startpoint_sd(:,2)<frameb.cs_start);
            end
            a=act_all(ind_event_in_this_fish,:);b=zeros(16,length(ind),length(ind_event_in_this_fish));
            for tt=1:length(ind)
                t=max(round(beh.startpoint_sd(ind(tt))/fs.behavior/fs.ca)-5,1):min(round(beh.startpoint_sd(ind(tt))/fs.behavior/fs.ca)+10,size(a,2));
                b(1:length(t),tt,:)=a(:,t)';
            end
            bb=squeeze(mean(b,2));
            plot(mean(bb,2),'color',color{kk},'linewidth',2);hold on;
            %shadedErrorBar([1:16],bb',{@(x) mean(x,1),@(x) std(x,[],1)},'lineprops',{color{kk}},'transparent',1,'patchSaturation',0.1); hold on
        end
        hold off;zz=zz+1;title([num2str(ii,'%02d')]);
    end
end
saveas(h,[outputpath '\motor trigger avg' ],'fig');
%% chose multi cluster num.
X=env_all.supervoxel;
for ii=unique(X(:,3))'
    id=find(X(:,3)==ii);
    X(id,3)=(ii-1)*8/0.66+1;
end
clrmap_iii=clrmap;cIX_iii=cIX;gIX_iii=gIX;a=[];
%num_clust=[1,5,15,18,19,20,22,25,30,32,33,35,38:41];num_clust=setdiff(unique(gIX_iii)',[1,5,15,18,19,20,22,25,30,32,33,35,38:41]);
num_clust=[4,10,20,21,23];%num_clust(num_clust==7)=[];
clrmap_iii=[(bl_fraction-bn_fraction)./(bl_fraction+bn_fraction)]'; clrmap_iii(find(abs(clrmap_iii)<0.2))=nan;
h_test_clus(num_clust);
table(h_acq_clus(num_clust,:));
for kk=[2,4,6,13];unique(index_all)';
    [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
    load(envbatch{kk});
    act_trace=act_all;
    ind_add=[];
    for ii=num_clust
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
       %ind_cls_in_this_fish( find(gIX_iii(i_cIX)~=ii))=[];
    end
    ind_cls_in_this_fish=i_cIX(ind_add);
    colorind=unique(gIX_iii(ind_cls_in_this_fish));
    
%     a=a+length(ind_cls_in_this_fish);
%     a=[a;cIX_iii(ind_cls_in_this_fish)];
    %h=pushbutton_popupplot_Callback(act_trace,cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),clrmap_iii(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
    h=DrawTiledPics_zyq_20190530(cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),[1:size(act_trace,1)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap_iii);
    %saveas(h,[savepath_all '\fish' num2str(kk,'%02d') '_' num2str(num_clust,'%02d')],'fig');
%     figure,        
%     %scatter3(X(ind_event_in_this_fish,1),X(ind_event_in_this_fish,2),X(ind_event_in_this_fish,3),8,[0.5 0.5 0.5],'filled');hold on;axis equal;
%     scatter3(X(cIX_iii(ind_cls_in_this_fish),1),X(cIX_iii(ind_cls_in_this_fish),2),X(cIX_iii(ind_cls_in_this_fish),3),12,clrmap_iii(gIX_iii(ind_cls_in_this_fish),:),'filled');axis equal;
%     axis equal;colorbar;colormap('hot');
%     set(gca,'color','k','visible','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[],'clim',[0.2 1])
%     view([-180 90]);
end
%% chose signle cluster num.
clrmap_iii=clrmap;cIX_iii=cIX;gIX_iii=gIX;a=[];zz=1;
figure,
num_clust=[4,10,20,21,23];clrmap_iii=1:length(num_clust);
num_clust=setdiff(unique(gIX_iii)',[1,5,15,18,19,20,22,25,30,32,33,35,38:41]);[11,15 17 21 22];
for ii=num_clust
    for kk=2;%unique(index_all)'
        [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
        load(envbatch{kk});
        %act_trace=act_all_acq';
        ind_cls_in_this_fish=i_cIX;ind_add=[];
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
        %ind_cls_in_this_fish( find(gIX_iii(i_cIX)~=ii))=[];
        ind_cls_in_this_fish=i_cIX(ind_add);
        colorind=unique(gIX_iii(ind_cls_in_this_fish));
        %a=a+length(ind_cls_in_this_fish);
        %a=[a;cIX_iii(ind_cls_in_this_fish)];
        %h=pushbutton_popupplot_Callback(act_trace,cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),clrmap_iii(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
        %h=DrawTiledPics_zyq_20190530(cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),[1:size(A,3)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap_iii);
        %saveas(h,[savepath_all '\fish' num2str(kk,'%02d') '_' num2str(ii,'%02d')],'fig');
        subplot(ceil(length(num_clust)/4),4,zz),xlim([1 env.height]);ylim([1 env.width]);
        scatter3(X(ind_event_in_this_fish,1),X(ind_event_in_this_fish,2),X(ind_event_in_this_fish,3),8,[0.5 0.5 0.5],'filled');hold on;axis equal;
        scatter3(X(cIX_iii(ind_cls_in_this_fish),1),X(cIX_iii(ind_cls_in_this_fish),2),X(cIX_iii(ind_cls_in_this_fish),3),12,clrmap_iii(gIX_iii(ind_cls_in_this_fish),:),'filled');hold on;
%         scatter(X(ind_event_in_this_fish,1),X(ind_event_in_this_fish,2),8,[0.5 0.5 0.5],'filled');hold on;axis equal;
%         scatter(X(cIX_iii(ind_cls_in_this_fish),1),X(cIX_iii(ind_cls_in_this_fish),2),12,clrmap_iii(gIX_iii(ind_cls_in_this_fish),:),'filled');hold on;
        axis equal;
        title(['Clust num ' num2str(ii,'%02d')]);
        view([-90 90]);
        zz=zz+1;
        %set(gca,'visible','off'); grid off;
    end
end
%% typical case
ind=cIX_iii(find(gIX_iii==4));ind=ind(202);
% figure,plot(squeeze(act_all_hab_mean(:,:,ind)));hold on;
% plot(squeeze(act_all_tst_mean(:,:,ind)));
figure,plot_rawtrace_trials(mean(act_all(ind,:),1)',[],fs,frame,trial,[],1);
c=corr(act_all(ind,:)',mean(act_all(ind,:),1)');figure,hist(c)
i_ind=find(c>0.8);i_ind=i_ind(100);
figure,plot_rawtrace_trials(mean(act_all(ind(i_ind),:),1)',[],fs,frame,trial,[],1);