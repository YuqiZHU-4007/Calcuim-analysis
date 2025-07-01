%stim regress(emerged population) & motion regress
%refer to Chen,2018,Neuron
% change GetFishFreq.m!!!!!!(related the shape of regressor)
%zyq
clc;clear all;close all;
set(0,'defaultfigurecolor','w');
%singel fish
%behaviorpath='E:\A_Data_lightsheet\Data_huc\20180509_huc\fish2\20180509_fish2_huc';
[actname,actpath]=uigetfile('','act_aftprocess.mat');
[envname,envpath]=uigetfile(actpath,'env.mat');
[bahname,behpath]=uigetfile(actpath,'behav');
[paraname,parapath]=uigetfile(actpath,'para.mat');
load([actpath actname]);
load([envpath envname]);
load([behpath bahname]);
load([parapath paraname]);
savepath=checkpath([actpath 'save\']);
%% batch
[n,p]=uigetfile('/*.txt','path');load([p,n]);
[n,p]=uigetfile('/*.txt','path');load([p,n]);
%load('G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\path.mat');
%load('K:\5.Calcium summary subset of G\huc_nonlearner_cutmove_20190921\path\path.mat')
%correct actpath
name='activities_aft_process';
for ii=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{ii});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{ii}=fullfile(fpath,[fname,ext]);
        savepathi{ii}=checkpath([fpath,'\save\']);
    end
end
savepath_eps=uigetdir(p);
%savepath_eps='G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\regress\';
%% motor regression:motor seeds
for batchi=1:length(path)
    load(actbatch{ii});
    load(behavbatch{ii});
    load(envbatch{ii});
    load(parabatch{ii});
    savepath=savepathi{ii};
    nCells_total = size(activities_preCS_dfdf_aftcorrect,2);prct_const = 2 ;topN = round(prct_const/100 * nCells_total);%para
    %分离spon and CS induce
    motor_CS=re_start_end_point(find(re_start_end_point(:,2,:)>= frameb.cs_start & re_start_end_point(:,2,:)< frameb.cs_end),:);
    motor_spon=re_start_end_point(find(re_start_end_point(:,2,:)< frameb.cs_start | re_start_end_point(:,2,:)>= frameb.cs_end),:);
    for hh=1:3
        switch hh
            case 1 %全长，所有shaking
                y=delta_c_bef1; figure,plot(y);
            case 2 %only CS
                y= zeros(size(delta_c_bef1));
                ind=[(motor_CS(:,1)-1)*frameb.per_cycle+motor_CS(:,2) (motor_CS(:,1)-1)*frameb.per_cycle+motor_CS(:,3)];
                for jj=1:size(ind,1)
                    y(ind(jj,1):ind(jj,2))=delta_c_bef1(ind(jj,1):ind(jj,2));
                end
                figure,plot(abs(y));hold on;
                a=[0 20];
                patch1=patch([[frameb.cs_start:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_start:frameb.per_cycle:size(y',1)]']',...
                    repmat([min(a) min(a) max(a) max(a)],length([frameb.cs_start:frameb.per_cycle:size(y',1)]),1)',...
                    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
            case 3 %only Spon
                y= zeros(size(delta_c_bef1));
                ind=[(motor_spon(:,1)-1)*frameb.per_cycle+motor_spon(:,2) (motor_spon(:,1)-1)*frameb.per_cycle+motor_spon(:,3)];
                for jj=1:size(ind,1)
                    y(ind(jj,1):ind(jj,2))=delta_c_bef1(ind(jj,1):ind(jj,2));
                end
                figure,plot(abs(y));hold on;
                a=[0 20];
                patch1=patch([[frameb.cs_start:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
                    [frameb.cs_start:frameb.per_cycle:size(y',1)]']',...
                    repmat([min(a) min(a) max(a) max(a)],length([frameb.cs_start:frameb.per_cycle:size(y',1)]),1)',...
                    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        end
        time_be=0:1/fs.behavior:(length(y)-1)/fs.behavior;
        time_ca=0:fs.ca:(length(activities_preCS_dfdf_aftcorrect(:,1))-1)*fs.ca;
        
        behavior=zeros(5,length(time_be));behavior(1,:)=abs(y);
        %behavior_down=zeros(5,length(time_ca));behavior_down(1,:)=interp1(time_be,behavior(1,:),time_ca);
        [regressors,~,~,~] = GetMotorRegressor(behavior,6);
        figure,plot(time_be,behavior(1,:));hold on;
        plot(time_be,regressors(1,1).im,'r')
        regressors_motor_down = interp1(time_be,regressors(1,1).im,time_ca);
        figure,plot(time_ca,regressors_motor_down);
        [~,mtRescorr] = MotorSourceCorrelation(activities_preCS_dfdf_aftcorrect',[],regressors_motor_down);
        figure,hist(mtRescorr)
        showspv= mapback(mtRescorr, env.supervoxel,[env.height env.width env.depth],1:length(mtRescorr));
        show_spv_GUI(showspv);
        
        [~,IX] = sort(mtRescorr,'descend');
        x0 = mtRescorr(IX(topN));x0=max(x0,0.55);
        
        ind_motor_output=find(mtRescorr> x0);
        motor_output=mean(activities_preCS_dfdf_aftcorrect(:,ind_motor_output)',1);
        figure,plot(motor_output);hold on;scatter(startpoint_sd/fs.behavior/fs.ca,-0.01*ones(size(startpoint_sd)),12,'k','^','filled');
        a=[-0.02 0.15];
        patch1=patch([[frame.cs_start:frame.per_cycle:size(motor_output',1)]'...
            [frame.cs_end:frame.per_cycle:size(motor_output',1)]'...
            [frame.cs_end:frame.per_cycle:size(motor_output',1)]'...
            [frame.cs_start:frame.per_cycle:size(motor_output',1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_start:frame.per_cycle:size(motor_output',1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        xlim([1 length(motor_output)])
        showspv  = mapback_show_zyq_20190729(ones(size(ind_motor_output')), env.supervoxel,ind_motor_output,env.vol,'dot');
        h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'motor seed\'],'tif');%saveas(h,[savepath,'motor seed.fig']);
        
        [~,mtRescorr] = MotorSourceCorrelation(activities_preCS_dfdf_aftcorrect',[],motor_output);
        figure,hist(mtRescorr);
        thr=mean(mtRescorr)+2*std(mtRescorr)
        line([thr,thr],[0 2000]);
        showspv= mapback(mtRescorr, env.supervoxel,[env.height env.width env.depth],1:length(mtRescorr));
        h=show_spv_GUI(showspv);seqwrite(showspv(:,:,1:24),[savepath 'motor corr mapping\'],'tif');%saveas(h,[savepath,'motor corr mapping.fig']);
        
        % ind_motor_output=find(mtRescorr>0.85);
        % motor_output=mean(activities_preCS_dfdf_aftcorrect(:,ind_motor_output)',1);
        % figure,plot(motor_output);
        % showspv  = mapback_show_zyq_20190729(ones(size(ind_motor_output')), env.supervoxel,ind_motor_output,env.vol,'dot');
        % h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'motor seed2\'],'tif');
        
        ind_pass_motor=find(mtRescorr>thr);
        showspv= mapback(mtRescorr(ind_pass_motor), env.supervoxel,[env.height env.width env.depth],ind_pass_motor);
        h=show_spv_GUI(showspv);seqwrite(showspv(:,:,1:24),[savepath 'motor corr mapping_pass_' num2str(thr) '\'],'tif');
        
        %% decomposition
        trial_ind=trial.test(2):trial.test(3);%trial.hab(2):trial.hab(3);%
        frameb_ind=frameb.cs_start:frameb.cs_end;
        frame_ind=frame.cs_start:frame.cs_end;
        stimCS=zeros(1,frame.per_cycle*length(trial_ind));
        for ii=1:size(stimCS,2)/frame.per_cycle
            stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=1;%23
        end
        motor_output_rep=reshape(motor_output,frame.per_cycle,[]);
        Stim_motor_trial=unique(re_startpoint_sd(find(re_startpoint_sd(:,1)>= trial_ind(1) & re_startpoint_sd(:,1)<= trial_ind(end) & re_startpoint_sd(:,2)>= frameb_ind(1) & re_startpoint_sd(:,2)<= frameb_ind(end)),1));
        Motor_trial=unique(re_startpoint_sd(find(re_startpoint_sd(:,1)>= trial_ind(1) & re_startpoint_sd(:,1)<= trial_ind(end) & (re_startpoint_sd(:,2)< frameb_ind(1) | re_startpoint_sd(:,2)> frameb_ind(end))),1));
        
        motor_tAvr=zeros(frame.per_cycle,length(trial_ind));
        motor_tRes=zeros(frame.per_cycle,length(trial_ind));
        for ii=1:length(trial_ind)
            if ~isempty(intersect(trial_ind(ii),Stim_motor_trial))
                motor_tAvr(:,ii)=[mean(motor_output_rep(1:frame_ind(1)-1,Stim_motor_trial),2);motor_output_rep(frame_ind(1):end,trial_ind(ii))];
            else
                motor_tAvr(:,ii)=mean(motor_output_rep(:,trial_ind),2);
                %          [mean(motor_output_rep(1:frame_ind(1)-1,setdiff(trial_ind,Stim_motor_trial)),2);...
                %              mean(motor_output_rep(frame_ind,Stim_motor_trial),2);...
                %              mean(motor_output_rep(frame_ind(end)+1:end,setdiff(trial_ind,Stim_motor_trial)),2)];
            end
            if ~isempty(intersect(trial_ind(ii),Motor_trial))
                motor_tRes(:,ii)=[motor_output_rep(1:frame_ind(1)-1,trial_ind(ii));mean(motor_output_rep(frame_ind,Motor_trial),2);motor_output_rep(frame_ind(end)+1:end,trial_ind(ii))];
            else
                motor_tRes(:,ii)=mean(motor_output_rep(:,trial_ind),2);
            end
        end
        motor_tAvr=reshape(motor_tAvr,1,[])';motor_tRes=reshape(motor_tRes,1,[])';
        
        %motor_tAvr=repmat(mean(motor_output_rep(:,Stim_motor_trial),2),length(trial_ind),1);
        %motor_output_rep=motor_output_rep(:,trial_ind);
        %motor_tRes=reshape(motor_output_rep,1,[])'-motor_tAvr;
        %behavior=reshape(motor_output_rep,1,[]);[~,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace_zyq_20190730(behavior,frame.per_cycle);
        
        figure,
        subplot(4,1,1),plot(stimCS);box off;ylim([0 2]);set(gca,'visible','off');%xlim([0 2300]);
        subplot(4,1,2),plot(reshape(motor_output_rep,1,[]));hold on;
        scatter(startpoint_sd/fs.behavior/fs.ca,-0.01*ones(size(startpoint_sd)),12,'k','^','filled');xlim([(trial_ind(1)-1)*frame.per_cycle+1 (trial_ind(end))*frame.per_cycle]);
        title('motor output');ylim([-0.05 0.15]);
        subplot(4,1,3),plot(motor_tAvr);title('motor avg');ylim([-0.05 0.15]);
        subplot(4,1,4),plot(motor_tRes);title('motor res');ylim([-0.05 0.15]);
        %% stim-independent and motor dependent
        M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
        reg_thres2=0.5;
        [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace_zyq_20190730(M,frame.per_cycle);
        b_stim = sqrt(abs(var(Data_tAvr')./var(Data_p'))); % b_stim
        mtAvrcorr = corr(motor_tAvr,Data_p');
        mtRescorr = corr(motor_tRes,Data_p');
        stimcorr = b_stim;
        [~,IX] = sort(mtRescorr,'descend');
        y0 = mtRescorr(IX(topN));y0=max(y0,reg_thres2);
        IX_passY = find(mtRescorr>=y0);
        IX_failY = find(mtRescorr<y0);
        [~,IX] = sort(mtAvrcorr,'descend');
        x0 = mtAvrcorr(IX(topN));x0=max(x0,reg_thres2);
        IX_passX = find(mtAvrcorr>=x0);
        IX_failX = find(mtAvrcorr<x0);
        IX_int = intersect(IX_passX,IX_passY);
        IX_Xonly = setdiff(IX_passX,IX_int);
        IX_Yonly = setdiff(IX_passY,IX_int);
        h = figure('Position',[500,500,300,400]); hold on
        %plot([0:0.1:1],[0:0.1:1],'--','color',[0.5 0.5 0.5]);
        scatter(mtAvrcorr,mtRescorr,3,[0.5,0.5,0.5],'filled');
        scatter(mtAvrcorr(IX_Yonly),mtRescorr(IX_Yonly),3,[0.9,0.2,0.9],'filled');
        scatter(mtAvrcorr(IX_Xonly),mtRescorr(IX_Xonly),3,[0.3,0.8,0],'filled');
        scatter(mtAvrcorr(IX_int),mtRescorr(IX_int),3,[0,0.3,0.9],'filled');
        xlabel('tAvr corr.');ylabel('tRes corr.');
        axis equal
        set(gca,'XTick',0:0.5:1,'YTick',-0.5:0.5:1);
        figure,
        subplot(4,1,1),plot(stimCS);box off;ylim([0 2]);set(gca,'visible','off');%xlim([0 2300]);
        subplot(4,1,2),plot(mean(M(IX_Yonly,:)));title('motor avg');hold on;
        %scatter(startpoint_sd/fs.behavior/fs.ca,-0.01*ones(size(startpoint_sd)),12,'k','^','filled');xlim([(trial_ind(1)-1)*frame.per_cycle+1 (trial_ind(end))*frame.per_cycle]);
        ylim([-0.05 0.15]);
        subplot(4,1,3),plot(mean(M(IX_Xonly,:)));title('motor res');ylim([-0.05 0.15]);
        subplot(4,1,4),plot(mean(M(IX_int,:)));title('overlaped');ylim([-0.05 0.15]);
        
        for ii=1
            %     switch ii
            %         case 1
            %             ind=IX_Yonly;clr=[0.9,0.2,0.9];
            %         case 2
            %             ind=IX_Xonly;clr=[0.3,0.8,0];
            %         case 3
            %             ind=IX_int;clr=[0,0.3,0.9];
            %     end
            ind=[IX_Yonly IX_Xonly IX_int];
            clr1=repmat([0.9,0.2,0.9],length(IX_Yonly),1);
            clr2=repmat([0.3,0.8,0],length(IX_Xonly),1);
            clr3=repmat([0,0.3,0.9],length(IX_int),1);
            clr=cat(1,clr1,clr2,clr3);clr=reshape(clr,length(ind),1,3);
            showspv=  mapback_show_zyq_20190729(clr, env.supervoxel,ind,env.vol(:,:,1:24),'dot');
            h=show_spv_GUI(showspv);
            seqwrite(showspv,[savepath 'motor decomposition corr mapping\'],'tif');
            %saveas(h,[savepath,'motor decomposition corr mapping.fig']);
        end
    end
end
%% stim regression
close all;
colorCS=[0.5 0.5 0.9];
reg_thres=0.01;%0.75
reg_thres2=0.5;%0.65
prct_const1 = 2 ;% top x%
prct_const2 = 50 ;% top x%
fishnum=2;
isautothr=false;
mulriple_thr_sd=3;
stAvrcorr_all_raw={};
for batchi=1:size(actbatch,1)
    load(actbatch{batchi});
    load(behavbatch{batchi});
    load(envbatch{batchi});
    load(parabatch{batchi});
    
    frame_ind_end=frame.us_start-1;%！！！！！！！！！！！！
    nCells_total = size(activities_preCS_dfdf_aftcorrect,2);%para
    trial_ind=1:trial.acq_block_trial;
    stimCS=zeros(1,frame.per_cycle*length(trial_ind));
    for ii=1:size(stimCS,2)/frame.per_cycle
        stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=1;%23
    end
    %regressors = GetStimRegressor(stim,fishset,i_fish);
    stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=stimCS;
    [stregressors,~,~,~] = GetMotorRegressor(stim,fishnum);
    %     figure,plot(stimCS);ylim([0 2]);hold on;
    %     plot(stregressors(1,1).im,'r')
    stAvrcorr_all=[];
    for ii=1:trial.trial.acq_block_num+2
        if ii>=2 && ii<=trial.trial.acq_block_num+1
            frame_ind=1:frame_ind_end;isus=true;
        else
            frame_ind=1:frame.per_cycle;isus=false;
        end
        switch ii
            case 1
                trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
            case {num2str([2:trial.trial.acq_block_num+1]')}
                trial_ind=trial.acq(2)+trial.acq_block_trial*(ii-2):trial.acq(2)+trial.acq_block_trial*(ii-1)-1;
            case trial.trial.acq_block_num+2
                trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
        end
        stim_output=[stregressors(1,1).im];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
        M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
        M=reshape(M,[],frame.per_cycle,length(trial_ind));M=M(:,frame_ind,:);M=reshape(M,[],length(frame_ind)*length(trial_ind),1);
        [stimcorr,~] = MotorSourceCorrelation(M,stim_output,[]);
        %figure,hist(stimcorr)
        showspv= mapback(stimcorr, env.supervoxel,[env.height env.width env.depth],1:length(stimcorr));
        %h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '\'],'tif');
        %set(h,'name',['stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
        %close(h);
        %saveas(h,[savepath,'stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '.fig']);
        [~,IX] = sort(stimcorr,'descend');
        topN = round(prct_const1/100 * nCells_total);
        x0 = stimcorr(IX(topN));
        ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
        %         if reg_thres>0
        %             ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
        %         else
        %             ind_stim_output = (find(stimcorr<max(reg_thres,x0)))';
        %         end
        if x0>reg_thres
            tag1='TOP2\';
        else
            tag1=[num2str(reg_thres),'\'];
        end
        %savepath=checkpath([savepathi{batchi}, tag1]);
        stim_output=mean(M(ind_stim_output,:),1);
        %wIX = stimcorr(ind_stim_output); % weight, i.e. corr.coeff
        %figure,plot(mean(M(ind_stim_output,:),1));
        
        %showspv  = mapback_show_zyq_20190729(stimcorr(ind_stim_output)', env.supervoxel,ind_stim_output,env.vol,'dot');
        %h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'stim corr mapping_stim seed' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '\'],'tif');
        %set(h,'name',['stim corr mapping_stim seed' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
        %close(h);
        %saveas(h,[savepath,'stim corr mapping.fig']);
        [stim_tAvr,stim_tRes,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
        stimCS2=reshape(stimCS,[],length(trial_ind));stimCS2=stimCS2(frame_ind,:);stimCS2=reshape(stimCS2,1,[]);
        %     figure;
        %     subplot(4,1,1,ax),plot(stimCS2);box off;ylim([0 2]);set(gca,'visible','off');xlim([1 length(stimCS2)]);
        %     subplot(4,1,2,ax),plot(reshape(stim_output,1,[]));hold on;title('stim output');xlim([1 length(stimCS2)]);
        %     subplot(4,1,3,ax),plot(stim_tAvr);title('stim avg');xlim([1 length(stimCS2)]);
        %     subplot(4,1,4,ax),plot(stim_tRes);title('stim res');xlim([1 length(stimCS2)]);
        
        stAvrcorr = corr(stim_tAvr',M'); 
        [~,IX] = sort(stAvrcorr,'descend');
        topN = round(prct_const2/100 * nCells_total);
        x0 = stAvrcorr(IX(topN));
        thr=max(x0,reg_thres2);
        IX_passX = find(stAvrcorr>=thr);figure,scatter(env.supervoxel(IX_passX,1),env.supervoxel(IX_passX,2))
        showspv  = mapback_show_zyq_20190729(stimcorr(IX_passX)', env.supervoxel,IX_passX,env.vol,'dot');
        
        %h=show_spv_GUI(showspv);%seqwrite(showspv,[savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)) '\'],'gif');
        %set(h,'name',['stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
        %close(h);
        if x0>reg_thres2
            savepath=checkpath([savepath, 'TOP2\']);
            filename = [savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)),'.gif']; % Specify the output file name
        else
            savepath=checkpath([savepath, num2str(reg_thres2) '\']);
            filename = [savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)),'.gif']; % Specify the output file name;('_thr_',regexprep(num2str(thr2),'[.]','$'),)
        end
        for idx = 1:size(showspv,4)
            [A,map] = rgb2ind(showspv(:,:,:,idx),256);
            if idx == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
            end
        end
        %saveas(h,[savepath,'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)),'.fig']);
    end
    h1=figure('name',['Num.',num2str(batchi),'_Stim_seed']);
    set(gcf,'position',[50,200,2500,150]);
    a=[-0.05 0.1];
    for ii=1:7
        if ii>=2 && ii<=trial.trial.acq_block_num+1
            frame_ind=1:frame_ind_end;isus=true;
        else
            frame_ind=1:frame.per_cycle;isus=false;
        end
        switch ii
            case 1
                trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
            case {num2str([2:trial.trial.acq_block_num+1]')}
                trial_ind=trial.acq(2)+trial.acq_block_trial*(ii-2):trial.acq(2)+trial.acq_block_trial*(ii-1)-1;
            case trial.trial.acq_block_num+2
                trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
        end
        stim_output=[stregressors(1,1).im];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
        M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
        M=reshape(M,[],frame.per_cycle,length(trial_ind));M=M(:,frame_ind,:);M=reshape(M,[],length(frame_ind)*length(trial_ind),1);
        [stimcorr,~] = MotorSourceCorrelation(M,stim_output,[]);
        [~,IX] = sort(stimcorr,'descend');
        topN = round(prct_const1/100 * nCells_total);
        x0 = stimcorr(IX(topN));
        ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
%         if reg_thres>0
%             ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
%         else
%             ind_stim_output = (find(stimcorr<max(reg_thres,x0)))';
%         end
        stim_output=mean(M(ind_stim_output,:),1);
        [stim_tAvr,stim_tRes,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
        stimCS2=reshape(stimCS,[],length(trial_ind));stimCS2=stimCS2(frame_ind,:);stimCS2=reshape(stimCS2,1,[]);
        subplot(1,7,ii);
        x=stimCS2';
        patch1=patch([[frame.cs_start:length(frame_ind):size(x,1)]'...
            [min(frame.cs_end,frame_ind(end)):length(frame_ind):size(x,1)]'...
            [min(frame.cs_end,frame_ind(end)):length(frame_ind):size(x,1)]'...
            [frame.cs_start:length(frame_ind):size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_start:length(frame_ind):size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        p1=plot(reshape(stim_output,1,[]),'k','linewidth',1.5);hold on;
        p2=plot(stim_tAvr,'r','linewidth',1);hold on;
        p3=plot(stim_tRes,'b','linewidth',1);hold on;
        box off;ylim(a);xlim([1 length(stimCS2)]);
        if ii==7,legend([patch1 p1,p2,p3],{'Stim.','stim output','stim tAvr','stim tRes'});end
    end
    saveas(h1,[savepath_eps 'Num.',num2str(batchi),'_Stim_seed' '.eps'],'epsc');

    h2=figure('name',['Num.',num2str(batchi),'_Ave.act & hist_of_StimAveCorr']);
    set(gcf,'position',[50,200,2500,330]);
    for ii=1:trial.acq_block_num+2
        if ii>=2 && ii<=trial.trial.acq_block_num+1
            frame_ind=1:frame_ind_end;isus=true;
        else
            frame_ind=1:frame.per_cycle;isus=false;
        end
        switch ii
            case 1
                trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
            case trial.acq_block_num+2
                trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
            otherwise
                trial_ind=trial.acq(2)+trial.acq_block_trial*(ii-2):trial.acq(2)+trial.acq_block_trial*(ii-1)-1;
        end
        stim_output=[stregressors(1,1).im];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
        M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
        M=reshape(M,[],frame.per_cycle,length(trial_ind));M=M(:,frame_ind,:);M=reshape(M,[],length(frame_ind)*length(trial_ind),1);
        [stimcorr,~] = MotorSourceCorrelation(M,stim_output,[]);
        
        [~,IX] = sort(stimcorr,'descend');
        topN = round(prct_const1/100 * nCells_total);
        x0 = stimcorr(IX(topN));
        ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
        %         if reg_thres>0
        %             ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
        %         else
        %             ind_stim_output = (find(stimcorr<max(reg_thres,x0)))';
        %         end
        stim_output=mean(M(ind_stim_output,:),1);
        
        [stim_tAvr,stim_tRes,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
        stimCS2=reshape(stimCS,[],length(trial_ind));stimCS2=stimCS2(frame_ind,:);stimCS2=reshape(stimCS2,1,[]);
        
        stAvrcorr = corr(stim_tAvr',M');
        stAvrcorr_all(ii,:)=stAvrcorr;%！！！
        [~,IX] = sort(stAvrcorr,'descend');
        topN = round(prct_const2/100 * nCells_total);
        x0 = stAvrcorr(IX(topN));
        thr=max(x0,reg_thres2);
        subplot(2,trial.acq_block_num+2,ii+trial.acq_block_num+2);hist(stAvrcorr);hold on;line([thr thr],[0 20000],'linewidth',2,'color','r');
        IX_passX = find(stAvrcorr>=thr);
        x=mean(activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,IX_passX),2);
        a=[-0.005 0.1];
        subplot(2,trial.acq_block_num+2,ii);
        patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_start:frame.per_cycle:size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        if isus
            l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
        end
        plot(x,'linewidth',1.5);hold on;
        ylim(a);
    end
    saveas(h2,[savepath_eps 'Num.',num2str(batchi),'_Ave.act & hist_of_StimAveCorr' '.eps'],'epsc');

    h3=figure('name',['Num.',num2str(batchi),'_hist_of_StimAveCorr']);
    c=[0 0 1;1 0 0];h=[];kk=1;
    for ii=[1,trial.acq_block_num+2]
        [~,IX] = sort(stAvrcorr_all(ii,:)','descend');
        topN = round(prct_const2/100 * nCells_total);
        x0 = stAvrcorr_all(ii,IX(topN));thr=max(x0,reg_thres2);
        h(kk,1)=histogram( stAvrcorr_all(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:));hold on
        line([thr thr],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
        kk=kk+1;
    end
    legend([h],{'Hab.','Tst.'});xlim([-1 1]);
    saveas(h3,[savepath_eps 'Num.',num2str(batchi),'_hist_of_StimAveCorr' '.eps'],'epsc');
    if isautothr
        auto_thres_bef_conditioing=mean(stAvrcorr_all(1,:))+mulriple_thr_sd*std(stAvrcorr_all(1,:));
        for ii=[1:7]
            if ii>=2 && ii<=6
                frame_ind=1:frame_ind_end;isus=true;
            else
                frame_ind=1:frame.per_cycle;isus=false;
            end
            switch ii
                case 1
                    trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
                case 2
                    trial_ind=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1;
                case 3
                    trial_ind=trial.acq(2)+trial.acq_block_trial*1:trial.acq(2)+trial.acq_block_trial*2-1;
                case 4
                    trial_ind=trial.acq(2)+trial.acq_block_trial*2:trial.acq(2)+trial.acq_block_trial*3-1;
                case 5
                    trial_ind=trial.acq(2)+trial.acq_block_trial*3:trial.acq(2)+trial.acq_block_trial*4-1;
                case 6
                    trial_ind=trial.acq(2)+trial.acq_block_trial*4:trial.acq(3);
                case 7
                    trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
            end
            stim_output=[stregressors(1,1).im];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
            M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
            M=reshape(M,[],frame.per_cycle,length(trial_ind));M=M(:,frame_ind,:);M=reshape(M,[],length(frame_ind)*length(trial_ind),1);
            [stimcorr,~] = MotorSourceCorrelation(M,stim_output,[]);
            %figure,hist(stimcorr)
            showspv= mapback(stimcorr, env.supervoxel,[env.height env.width env.depth],1:length(stimcorr));
            %h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '\'],'tif');
            %set(h,'name',['stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
            %close(h);
            %saveas(h,[savepath,'stim corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '.fig']);
            [~,IX] = sort(stimcorr,'descend');
            topN = round(prct_const1/100 * nCells_total);
            x0 = stimcorr(IX(topN));
            ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
            %         if reg_thres>0
            %             ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
            %         else
            %             ind_stim_output = (find(stimcorr<max(reg_thres,x0)))';
            %         end
            if x0>reg_thres
                tag1='TOP2\';
            else
                tag1=[num2str(reg_thres),'\'];
            end
            %savepath=checkpath([savepathi{batchi}, tag1]);
            stim_output=mean(M(ind_stim_output,:),1);
            %wIX = stimcorr(ind_stim_output); % weight, i.e. corr.coeff
            %figure,plot(mean(M(ind_stim_output,:),1));
            
            showspv  = mapback_show_zyq_20190729(stimcorr(ind_stim_output)', env.supervoxel,ind_stim_output,env.vol,'dot');
            %h=show_spv_GUI(showspv);seqwrite(showspv,[savepath 'stim corr mapping_stim seed' num2str(trial_ind(1)) '-' num2str(trial_ind(end)), '\'],'tif');
            %set(h,'name',['stim corr mapping_stim seed' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
            %close(h);
            %saveas(h,[savepath,'stim corr mapping.fig']);
            [stim_tAvr,stim_tRes,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
            stimCS2=reshape(stimCS,[],length(trial_ind));stimCS2=stimCS2(frame_ind,:);stimCS2=reshape(stimCS2,1,[]);
            %     figure;
            %     subplot(4,1,1,ax),plot(stimCS2);box off;ylim([0 2]);set(gca,'visible','off');xlim([1 length(stimCS2)]);
            %     subplot(4,1,2,ax),plot(reshape(stim_output,1,[]));hold on;title('stim output');xlim([1 length(stimCS2)]);
            %     subplot(4,1,3,ax),plot(stim_tAvr);title('stim avg');xlim([1 length(stimCS2)]);
            %     subplot(4,1,4,ax),plot(stim_tRes);title('stim res');xlim([1 length(stimCS2)]);
            
            stAvrcorr = corr(stim_tAvr',M'); 
            [~,IX] = sort(stAvrcorr,'descend');
            topN = round(prct_const2/100 * nCells_total);
            x0 = stimcorr(IX(topN));
            IX_passX = find(stAvrcorr>=auto_thres_bef_conditioing);
            showspv  = mapback_show_zyq_20190729(stimcorr(IX_passX)', env.supervoxel,IX_passX,env.vol,'dot');
            %h=show_spv_GUI(showspv);%seqwrite(showspv,[savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)) '\'],'gif');
            %set(h,'name',['stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end))]);
            %close(h);
            savepath=checkpath([savepath, 'autothr',num2str(auto_thres_bef_conditioing),'\']);
            filename = [savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)),'.gif']; % Specify the output file name
            for idx = 1:size(showspv,4)
                [A,map] = rgb2ind(showspv(:,:,:,idx),256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
                end
            end
            %saveas(h,[savepath,'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)),'.fig']);
        end
        
        h4=figure('name',['Num.',num2str(batchi),'_Ave.act & hist_of_StimAveCorr_autothr']);
        set(gcf,'position',[50,200,2500,330]);
        for ii=1:7
            if ii>=2 && ii<=6
                frame_ind=1:frame_ind_end;isus=true;
            else
                frame_ind=1:frame.per_cycle;isus=false;
            end
            switch ii
                case 1
                    trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
                case 2
                    trial_ind=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1;
                case 3
                    trial_ind=trial.acq(2)+trial.acq_block_trial*1:trial.acq(2)+trial.acq_block_trial*2-1;
                case 4
                    trial_ind=trial.acq(2)+trial.acq_block_trial*2:trial.acq(2)+trial.acq_block_trial*3-1;
                case 5
                    trial_ind=trial.acq(2)+trial.acq_block_trial*3:trial.acq(2)+trial.acq_block_trial*4-1;
                case 6
                    trial_ind=trial.acq(2)+trial.acq_block_trial*4:trial.acq(3);
                case 7
                    trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
            end
            stim_output=[stregressors(1,1).im];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
            M=activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,:)';
            M=reshape(M,[],frame.per_cycle,length(trial_ind));M=M(:,frame_ind,:);M=reshape(M,[],length(frame_ind)*length(trial_ind),1);
            [stimcorr,~] = MotorSourceCorrelation(M,stim_output,[]);
            
            [~,IX] = sort(stimcorr,'descend');
            topN = round(prct_const1/100 * nCells_total);
            x0 = stimcorr(IX(topN));
            ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
            %         if reg_thres>0
            %             ind_stim_output = (find(stimcorr>max(reg_thres,x0)))';
            %         else
            %             ind_stim_output = (find(stimcorr<max(reg_thres,x0)))';
            %         end
            stim_output=mean(M(ind_stim_output,:),1);
            
            [stim_tAvr,stim_tRes,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
            stimCS2=reshape(stimCS,[],length(trial_ind));stimCS2=stimCS2(frame_ind,:);stimCS2=reshape(stimCS2,1,[]);
            
            stAvrcorr = corr(stim_tAvr',M');
            [~,IX] = sort(stAvrcorr,'descend');
            topN = round(prct_const2/100 * nCells_total);
            x0 = stAvrcorr(IX(topN));
            thr=auto_thres_bef_conditioing;
            subplot(2,7,ii+7);hist(stAvrcorr);hold on;line([thr thr],[0 20000],'linewidth',2,'color','r');
            IX_passX = find(stAvrcorr>=thr);
            x=mean(activities_preCS_dfdf_aftcorrect((trial_ind(1)-1)*frame.per_cycle+1:(trial_ind(end))*frame.per_cycle,IX_passX),2);
            a=[-0.005 0.1];
            subplot(2,7,ii);
            patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
                [frame.cs_end:frame.per_cycle:size(x,1)]'...
                [frame.cs_end:frame.per_cycle:size(x,1)]'...
                [frame.cs_start:frame.per_cycle:size(x,1)]']',...
                repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
                colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
            if isus
                l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
            end
            plot(x,'linewidth',1.5);hold on;
            ylim(a);
        end
        saveas(h4,[savepath_eps 'Num.',num2str(batchi),'_Ave.act & hist_of_StimAveCorr_autothr' '.eps'],'epsc');
        
        h5=figure('name',['Num.',num2str(batchi),'_hist_of_StimAveCorr_autothr']);
        c=GetColormap('hsv_new',5);%[0 0 1;1 0 0,];
        kk=1;h=[];
        for ii=[2:6]
            [~,IX] = sort(stAvrcorr_all(ii,:)','descend');
            h(kk,1)=histogram( stAvrcorr_all(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:));hold on
            line([auto_thres_bef_conditioing auto_thres_bef_conditioing],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
            kk=kk+1;
        end
        %legend([h],{'Acq.1','Acq.2','Acq.3','Acq.4','Acq.5'});
        legend([h],{'Hab.','Acq.1','Acq.2','Acq.3','Acq.4','Acq.5','Tst.'});
        %legend([h],{'Hab.','Tst.'});
        xlim([-1 1]);
        saveas(h5,[savepath_eps 'Num.',num2str(batchi),'_hist_of_StimAveCorr_autothr' '.eps'],'epsc');
    end
    stAvrcorr_all_raw{batchi,1}=stAvrcorr_all;
end
save([savepath_eps 'stAvrcorr_all_raw.mat'],'stAvrcorr_all_raw','reg_thres','reg_thres2','prct_const1','prct_const2','fishnum','isautothr','mulriple_thr_sd');
close all;
% figure,
% colormap hot;
% imagesc(stAvrcorr_all',[0 1]);colorbar;
% eucD = pdist(stAvrcorr_all','euclidean');
% clustTreeEuc = linkage(eucD,'average');
% cophenet(clustTreeEuc,eucD)
% [h,nodes] = dendrogram(clustTreeEuc,200);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = [];
[cidx,~] = kmeans(stAvrcorr_all',4,'replicates',5,'dist','sqeuclidean');
%[silh2,h] = silhouette(stAvrcorr_all',cidx,'sqeuclidean');
[~,scidx]=sort(cidx);
figure,colormap hot;
imagesc(stAvrcorr_all(:,scidx)',[0 1]);colorbar;

figure,hist(stAvrcorr_all(:))
stAvrcorr_all_b=stAvrcorr_all >= 0.75;
stAvrcorr_all_b(:,find(sum(stAvrcorr_all_b)==0))=[];
figure,colormap hot;
imagesc(stAvrcorr_all_b',[0 1]);colorbar;
% figure,colormap hot;
% imagesc(stAvrcorr_all_b(:,scidx)',[0 1]);colorbar;

% %sort
S=[];%score
w=repmat([1:size(stAvrcorr_all_b,1)],size(stAvrcorr_all_b,2),1)';
S=w.*stAvrcorr_all_b;
S=sum(S,1);
[~,I]=sort(S,'descend');
figure,
colormap hot;
imagesc(stAvrcorr_all_b(:,I)',[0 1]);colorbar;

R = corrcoef(stAvrcorr_all);
figure,imagesc(R);

figure,
plot([-1:0.5:1],[-1:0.5:1],'--','color',[0.5 0.5 0.5],'linewidth',2); hold on;
scatplot(stAvrcorr_all(1,:),stAvrcorr_all(2,:),[],[],[],[],2,[]);
axis equal;xlim([-1 1]);ylim([-1 1]);

%% shaking profile and detection
%每个trial最后一个block和后一个trial第一个block
frameb_ind=[frame.cs_start,frame.us_start-1];%%%%%%%%%%%%%%%%%%%!!!!!!!
trace=struct;
ind=(trial.hab(3)-1)*frame.per_cycle+frameb_ind(1):(trial.hab(3)-1)*frame.per_cycle+frameb_ind(2);
trace.hab=activities_preCS_dfdf_aftcorrect(ind,:);
ind=(trial.test(2)-1)*frame.per_cycle+frameb_ind(1):(trial.test(2)-1)*frame.per_cycle+frameb_ind(2);
trace.test=activities_preCS_dfdf_aftcorrect(ind,:);
indd=trial.acq(2):trial.acq_block_trial:trial.acq(3);ind=[];
for ii=1:length(indd)
    ind=[ind (indd(ii)-1)*frame.per_cycle+frameb_ind(1):(indd(ii)-1)*frame.per_cycle+frameb_ind(2)];
end
trace.acq=activities_preCS_dfdf_aftcorrect(ind,:);

%%shaking profile
CS_mov=re_startpoint(find(re_startpoint(:,2)>frameb.cs_start & re_startpoint(:,2)<frameb.cs_end & re_startpoint(:,1)>trial.acq(3)),:);
%CS_mov=re_startpoint(find(re_startpoint(:,2)<frameb.cs_start),:);
CS_mov(:,2)=CS_mov(:,2)/fs.behavior/fs.ca;
bin=10;act=[];
for ii=1:size(CS_mov,1)
    ind=CS_mov(ii,:);
    act(:,:,ii)=activities_preCS_dfdf_aftcorrect(ind(1)*frame.per_cycle+floor(ind(2))-bin:...
        min(ind(1)*frame.per_cycle+floor(ind(2))+bin,(ind(1)+1)*frame.per_cycle),:);
end
m_act=mean(act,3);
figure,plot(m_act);hold on;
line([bin bin],[-0.1 0.5],'color','r','linewidth',1.5,'linestyle','--');hold on;

for ii=1:size(m_act,2)
    a=m_act(:,:,ii);
    event(1,ii)=isevent_20190508(a,ref_win,area_win_hab,trial);
end
