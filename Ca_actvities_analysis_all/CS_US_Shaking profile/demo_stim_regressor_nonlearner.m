clc;clear all;close all;
set(0,'defaultfigurecolor','w');
load('E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_nonlearner_cutmove_20190921\path\path.mat')
name='activities_aft_process';
for ii=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{ii});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{ii}=fullfile(fpath,[fname,ext]);
        savepathi{ii}=checkpath([fpath,'\save\']);
    end
end
savepath_eps='K:\3.Poster_fig\part3\non-learner\';

close all;
colorCS=[0.5 0.5 0.9];
reg_thres=0.01;%0.75
reg_thres2=0.5;%0.65
prct_const1 = 2 ;% top x%
prct_const2 = 50 ;% top x%
fishnum=1;
isautothr=true;
mulriple_thr_sd=3;
stAvrcorr_all_raw={};
for batchi=1:size(actbatch,1)
    load(actbatch{batchi});
    load(behavbatch{batchi});
    load(envbatch{batchi});
    load(parabatch{batchi});
    
    frame_ind_end=frame.us_start;%ŁĄŁĄŁĄŁĄŁĄŁĄŁĄŁĄŁĄŁĄŁĄŁĄ
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
    for ii=1:7
        if ii>=2 || ii<=6
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
        savepath=checkpath([savepathi{batchi}, tag1]);
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
        IX_passX = find(stAvrcorr>=thr);
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
        if ii>=2 || ii<=6
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
    for ii=1:7
        if ii>=2 || ii<=6
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
        stAvrcorr_all(ii,:)=stAvrcorr;%ŁĄŁĄŁĄ
        [~,IX] = sort(stAvrcorr,'descend');
        topN = round(prct_const2/100 * nCells_total);
        x0 = stAvrcorr(IX(topN));
        thr=max(x0,reg_thres2);
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
    saveas(h2,[savepath_eps 'Num.',num2str(batchi),'_Ave.act & hist_of_StimAveCorr' '.eps'],'epsc');

    h3=figure('name',['Num.',num2str(batchi),'_hist_of_StimAveCorr']);
    c=[0 0 1;1 0 0];h=[];kk=1;
    for ii=[1,7]
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
            if ii>=2 || ii<=6
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
            savepath=checkpath([savepathi{batchi}, tag1]);
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
            if ii>=2 || ii<=6
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
        c=GetColormap('hsv_new',7);%[0 0 1;1 0 0,];
        kk=1;h=[];
        for ii=[1:7]
            [~,IX] = sort(stAvrcorr_all(ii,:)','descend');
            h(kk,1)=histogram( stAvrcorr_all(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:));hold on
            line([auto_thres_bef_conditioing auto_thres_bef_conditioing],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
            kk=kk+1;
        end
        legend([h],{'Hab.','Acq.1','Acq.2','Acq.3','Acq.4','Acq.5','Tst.'});
        %legend([h],{'Hab.','Tst.'});
        xlim([-1 1]);
        saveas(h5,[savepath_eps 'Num.',num2str(batchi),'_hist_of_StimAveCorr_autothr' '.eps'],'epsc');
    end
    stAvrcorr_all_raw{batchi,1}=stAvrcorr_all;
end
save([savepath_eps 'stAvrcorr_all_raw.mat'],'stAvrcorr_all_raw','reg_thres','reg_thres2','prct_const1','prct_const2','fishnum','isautothr','mulriple_thr_sd');
