clc;
clear all;
%% load all 
iscutmov=1;fraction_type=1;
savepath='X:\calcium data 20230224\';%'X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
%warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load(fullfile(Path{1}{1},'para.mat'));

pre_cs_time=[4 1];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca);
onset_win=frame.cs_start-4:frame.per_cycle;
A=struct;
path_all={};kk=1;
NUM_CR_UR_all=[];NUM_CR_UR_align_all=[];fraction_CR_UR_all=[];fraction_CR_UR_align_all=[];
CS_response_amp_all=[];CS_response_amp_align_all=[];region_CS_response_amp_all=[];region_CS_response_amp_align_all=[];
supervol_type_all={};supervol_type_align_all={};region_supervol_all={};region_supervol_align_all={};
supervoltemp_type_all={};supervoltemp_type_align_all={};region_supervoltemp_all={};region_supervoltemp_align_all={};
Response_amp_all=[];Response_amp_region_all=[];Response_amp_align_all=[];Response_amp_region_align_all=[];Response_amp_all_session_align_all=[];
Response_amp_region_smooth_all=nan(size(Response_amp_region_all));Response_amp_region_smooth_align_all=[];
start_index_all=[];start_index_align_all=[];
% A=setfield(A,'NUM_CR_UR',[]);A=setfield(A,'fraction_CR_UR',[]);
% A=setfield(A,'NUM_CR_UR_align',[]);A=setfield(A,'fraction_CR_UR_align',[]);
% A=setfield(A,'CS_response_amp',[]);A=setfield(A,'CS_response_amp_align',[]);
% A=setfield(A,'region_CS_response_amp',[]);A=setfield(A,'region_CS_response_amp_align',[]);
% A=setfield(A,'supervol_type',[]);A=setfield(A,'supervol_type_align',[]);A=setfield(A,'region_supervol',[]);A=setfield(A,'region_supervol_align',[]);
% A=setfield(A,'supervoltemp_type',[]);A=setfield(A,'supervoltemp_type_align',[]);A=setfield(A,'region_supervoltemp',[]);A=setfield(A,'region_supervoltemp_align',[]);
for groupi=1:4
    path=Path{groupi};
    for fishi=1:length(path)
        p=path{fishi};path_all{kk,1}=p;kk=kk+1;
        load(fullfile(p,'singletrial','CR_ind_summary_singletrial.mat'),'A_r','ind_all_CSUS_RESPONSIVE','Fraction_in_region_type','US_activation_response_cuttrail', 'CSon_activation_response','Index_in_region');
        load(fullfile(p, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo');
        load(fullfile(p, '/env.mat'));supervol=env.supervoxel(:,1:3);
        nn=[p(end-14:end-7),p(end-5:end-1)];
        supervol_temp=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        NUM_CR_UR=[];fraction_CR_UR=[];NUM_CR_UR_align=[];fraction_CR_UR_align=[];
        Response_amp=[];Response_amp_allsession=[];Response_amp_region=nan(frame.per_cycle,trial.total,30,size(ind_all_CSUS_RESPONSIVE,1));
        Response_amp_align=[];Response_amp_region_align=[]; Response_amp_align_allsession=[]; Response_amp_align_allsession1=[];
        CS_response_amp=[];CS_response_amp_align=[];region_CS_response_amp=[];region_CS_response_amp_align=[];
        supervol_type={};supervoltemp_type={};region_supervol={};region_supervoltemp={};
        supervol_type_align={};supervoltemp_type_align={};region_supervol_align={};region_supervoltemp_align={};
        start_index=nan(3,30,trial.total,size(ind_all_CSUS_RESPONSIVE,1));Response_amp_region_smooth=[];start_index_align=[];Response_amp_region_smooth_align=[];
        for typei=1:size(ind_all_CSUS_RESPONSIVE,1)
            
            for sessioni=1:size(Fraction_in_region_type{typei,iscutmov,fraction_type},2)
                
                NUM_CR_UR(sessioni,typei)=length(ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov})./size(supervol,1);
                fraction_CR_UR(sessioni,:,typei)=Fraction_in_region_type{typei,iscutmov,fraction_type}(:,sessioni);
                CS_response_amp(sessioni,typei)=squeeze(mean( mean(CSon_activation_response{sessioni,iscutmov}(:,:,ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov}),1,'omitnan'),3,'omitnan'));
                supervol_type{typei,sessioni}=supervol(ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov},:);
                supervoltemp_type{typei,sessioni}=supervol_temp(ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov},:);
                Response_amp(:,sessioni,typei)=squeeze(mean(A_r(:,sessioni,ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov}),3,'omitnan'));
                Response_amp_allsession(:,:,sessioni,typei)=squeeze(mean(A_r(:,:,ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov}),3,'omitnan'));
                for regioni=1:size(Fraction_in_region_type{typei,iscutmov,fraction_type},1)
                    a= Index_in_region{typei,sessioni,iscutmov,fraction_type}{regioni,:};
                    region_CS_response_amp(sessioni,regioni,typei)=squeeze(mean( mean(CSon_activation_response{sessioni,iscutmov}(:,:,a),1,'omitnan'),3,'omitnan'));
                    region_supervol{typei,regioni,sessioni}=supervol(a,:);
                    region_supervoltemp{typei,regioni,sessioni}=supervol_temp(a,:);
                    Response_amp_region(:,sessioni,regioni,typei)=squeeze(mean(A_r(:,sessioni,a),3,'omitnan'));
                end  
%                 figure, scatter3( supervol(:,1), supervol(:,2), supervol(:,3));hold on;
%                 scatter3(  region_supervol{typei,regioni,sessioni}(:,1), region_supervol{typei,regioni,sessioni}(:,2), region_supervol{typei,regioni,sessioni}(:,3))
                 a=squeeze(Response_amp_region(:,sessioni,:,typei));
                [ind,y,m,sd]=findonset(a,base_win,onset_win);
                start_index(:,:,sessioni,typei)=cat(1,ind,m,sd);
                Response_amp_region_smooth(:,sessioni,:,typei)=y;
            end
            for align_typei=1:3
                for session_align=1:size(align_win,2)
                    switch align_typei
                        case 1
                            a=align_win(:,session_align);
                        case 2
                            a=align_win_go(:,session_align);
                        case 3
                            a=align_win_nogo(:,session_align);
                    end
                    a(isnan(a))=[];
                    NUM_CR_UR_align(session_align,typei,align_typei)=mean(NUM_CR_UR(a,typei),1,'omitnan');
                    fraction_CR_UR_align(session_align,:,typei,align_typei)=mean(fraction_CR_UR(a,:,typei),1,'omitnan');
                    CS_response_amp_align(session_align,typei,align_typei)=mean(CS_response_amp(a,typei),1,'omitnan');
                    region_CS_response_amp_align(session_align,:,typei,align_typei)=mean(region_CS_response_amp(a,:,typei),1,'omitnan');
                    supervol_type_align{session_align,typei,align_typei}={supervol_type{typei,a}};
                    supervoltemp_type_align{session_align,typei,align_typei}={supervoltemp_type{typei,a}};
                    Response_amp_align(:,session_align,typei,align_typei)=mean(Response_amp(:,a,typei),2,'omitnan');
                    Response_amp_align_allsession1(:,:,session_align,typei,align_typei)=mean(Response_amp_allsession(:,:,a,typei),3,'omitnan');
                    
                    for regioni=1:size(region_supervol,2)
                        region_supervol_align{typei,regioni,session_align,align_typei}={region_supervol{typei,regioni,a}};
                        region_supervoltemp_align{typei,regioni,session_align,align_typei}={region_supervoltemp{typei,regioni,a}};
                        Response_amp_region_align(:,session_align,regioni,typei,align_typei)=mean(Response_amp_region(:,a,regioni,typei),2,'omitnan');
                        Response_amp_region_smooth_align(:,session_align,regioni,typei,align_typei)=mean(Response_amp_region_smooth(:,a,regioni,typei),2,'omitnan');
                    end
                    start_index_align(:,:,session_align,typei,align_typei)=mean(start_index(:,:,a,typei),3,'omitnan');
                end
                  for session_align=1:size(align_win,2)
                    switch align_typei
                        case 1
                            a=align_win(:,session_align);
                        case 2
                            a=align_win_go(:,session_align);
                        case 3
                            a=align_win_nogo(:,session_align);
                    end
                    a(isnan(a))=[];
                    Response_amp_align_allsession(:,session_align,:,typei,align_typei)=mean(Response_amp_align_allsession1(:,a,:,typei,align_typei),2,'omitnan');
                  end
            end
        end
        
        %
        a=NUM_CR_UR_all;a=cat(3,a,NUM_CR_UR);NUM_CR_UR_all=a;
        a=NUM_CR_UR_align_all;a=cat(4,a,NUM_CR_UR_align);NUM_CR_UR_align_all=a;
        a=fraction_CR_UR_all;a=cat(4,a,fraction_CR_UR);fraction_CR_UR_all=a;
        a=fraction_CR_UR_align_all;a=cat(5,a,fraction_CR_UR_align);fraction_CR_UR_align_all=a;
        a=CS_response_amp_all;a=cat(3,a,CS_response_amp);CS_response_amp_all=a;
        a=CS_response_amp_align_all;a=cat(4,a,CS_response_amp_align);CS_response_amp_align_all=a;
        a=region_CS_response_amp_all;a=cat(4,a,region_CS_response_amp);region_CS_response_amp_all=a;
        a=region_CS_response_amp_align_all;a=cat(5,a,region_CS_response_amp_align);region_CS_response_amp_align_all=a;
        a=supervol_type_all;a=cat(3,a,supervol_type);supervol_type_all=a;
        a=supervol_type_align_all;a=cat(4,a,supervol_type_align);supervol_type_align_all=a;
        a=region_supervol_all;a=cat(4,a,region_supervol);region_supervol_all=a;
        a=region_supervol_align_all;a=cat(5,a,region_supervol_align);region_supervol_align_all=a;
        a=supervoltemp_type_all;a=cat(3,a,supervoltemp_type);supervoltemp_type_all=a;
        a=supervoltemp_type_align_all;a=cat(4,a,supervoltemp_type_align);supervoltemp_type_align_all=a;
        a=region_supervoltemp_all;a=cat(4,a,region_supervoltemp);region_supervoltemp_all=a;
        a=region_supervoltemp_align_all;a=cat(5,a,region_supervoltemp_align);region_supervoltemp_align_all=a;
        a=Response_amp_all;a=cat(4,a,Response_amp);Response_amp_all=a;
        a=Response_amp_region_all;a=cat(5,a,Response_amp_region);Response_amp_region_all=a;
        a=Response_amp_align_all;a=cat(5,a,Response_amp_align);Response_amp_align_all=a;
        a=Response_amp_all_session_align_all;a=cat(6,a,Response_amp_align_allsession);Response_amp_all_session_align_all=a;
        a=Response_amp_region_align_all;a=cat(6,a,Response_amp_region_align);Response_amp_region_align_all=a;
        a=start_index_all;a=cat(5,a,start_index);start_index_all=a;
        a=Response_amp_region_smooth_all;a=cat(5,a,Response_amp_region_smooth);Response_amp_region_smooth_all=a;
        a=start_index_align_all;a=cat(6,a,start_index_align);start_index_align_all=a;
        a=Response_amp_region_smooth_align_all;a=cat(6,a,Response_amp_region_smooth_align);Response_amp_region_smooth_align_all=a;
        %
%         a=getfield(A,'NUM_CR_UR');a=cat(3,a,NUM_CR_UR);A=setfield(A,'NUM_CR_UR',a);
%         a=getfield(A,'NUM_CR_UR_align');a=cat(4,a,NUM_CR_UR_align);A=setfield(A,'NUM_CR_UR_align',a);
%         a=getfield(A,'fraction_CR_UR');a=cat(4,a,fraction_CR_UR);A=setfield(A,'fraction_CR_UR',a);
%         a=getfield(A,'fraction_CR_UR_align');a=cat(5,a,fraction_CR_UR_align);A=setfield(A,'fraction_CR_UR_align',a);
%         a=getfield(A,'CS_response_amp');a=cat(3,a,CS_response_amp);A=setfield(A,'CS_response_amp',a);
%         a=getfield(A,'CS_response_amp_align');a=cat(4,a,CS_response_amp_align);A=setfield(A,'CS_response_amp_align',a);
%         a=getfield(A,'region_CS_response_amp');a=cat(4,a,region_CS_response_amp);A=setfield(A,'region_CS_response_amp',a);
%         a=getfield(A,'region_CS_response_amp_align');a=cat(5,a,region_CS_response_amp_align);A=setfield(A,'region_CS_response_amp_align',a);
%         a=getfield(A,'supervol_type');a=cat(3,a,supervol_type);A=setfield(A,'supervol_type',a);
%         a=getfield(A,'supervol_type_align');a=cat(4,a,supervol_type_align);A=setfield(A,'supervol_type_align',a);
%         a=getfield(A,'region_supervol');a=cat(4,a,region_supervol);A=setfield(A,'region_supervol',a);
%         a=getfield(A,'region_supervol_align');a=cat(5,a,region_supervol_align);A=setfield(A,'region_supervol_align',a);
%         a=getfield(A,'supervoltemp_type');a=cat(3,a,supervoltemp_type);A=setfield(A,'supervoltemp_type',a);
%         a=getfield(A,'supervoltemp_type_align');a=cat(4,a,supervoltemp_type_align);A=setfield(A,'supervoltemp_type_align',a);
%         a=getfield(A,'region_supervoltemp');a=cat(4,a,region_supervoltemp);A=setfield(A,'region_supervoltemp',a);
%         a=getfield(A,'region_supervoltemp_align');a=cat(5,a,region_supervoltemp_align);A=setfield(A,'region_supervoltemp_align',a);
    end 
end
% save([savepath,'/Summary_all_CRUR_singletrial.mat'],'A','path_all','Path','iscutmov','fraction_type','-v7.3');
save([savepath,'/Summary_all_CRUR_singletrial.mat'],'Path','path_all','iscutmov','fraction_type','base_win','onset_win',...
    'NUM_CR_UR_all','NUM_CR_UR_align_all','fraction_CR_UR_all','fraction_CR_UR_align_all',...
'CS_response_amp_all','CS_response_amp_align_all','region_CS_response_amp_all','region_CS_response_amp_align_all',...
'supervol_type_all','supervol_type_align_all','region_supervol_all','region_supervol_align_all',...
'supervoltemp_type_all','supervoltemp_type_align_all','region_supervoltemp_all','region_supervoltemp_align_all',...
'Response_amp_all','Response_amp_region_all','Response_amp_align_all','Response_amp_region_align_all',...
'Response_amp_region_smooth_all','Response_amp_region_smooth_align_all','start_index_all','start_index_align_all','-v7.3');

%%
clc;
clear all;
iscutmov=1;fraction_type=1;
savepath='X:\calcium data 20230224\';%
%savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
%warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load(fullfile(Path{1}{1},'para.mat'));
Response_amp_all_session_align_all=[];Response_sd_all_session_align_all=[];
Response_amp_region_all_session_align_all=[];Response_sd_region_all_session_align_all=[];
for groupi=1
    path=Path{groupi};
    for fishi=1:length(path)
        p=path{fishi};
        load(fullfile(p,'singletrial','CR_ind_summary_singletrial.mat'),'A_r','ind_all_CSUS_RESPONSIVE','Fraction_in_region_type','Index_in_region');
        load(fullfile(p, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo');
        load(fullfile(p, '/env.mat'));supervol=env.supervoxel(:,1:3);
        nn=[p(end-14:end-7),p(end-5:end-1)];
        Response_amp_allsession=[];Response_amp_align_allsession=[]; Response_amp_align_allsession1=[];Response_sd_align_allsession1=[];
        Response_sd_allsession=[];Response_sd_align_allsession=[];
        Response_amp_region_allsession=[]; Response_amp_region_allsession1=[];Response_amp_region_align_allsession=[];
        Response_sd_region_allsession=[]; Response_sd_region_allsession1=[];Response_sd_region_align_allsession=[];
        for typei=[1,5,7,8,16,17]
            for sessioni=1:size(Fraction_in_region_type{typei,iscutmov,fraction_type},2)
                Response_amp_allsession(:,:,sessioni,typei)=squeeze(mean(A_r(:,:,ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov,1}),3,'omitnan'));
                Response_sd_allsession(:,:,sessioni,typei)=squeeze(std(A_r(:,:,ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov,1}),[],3,'omitnan'));
                for regioni=1:size(Fraction_in_region_type{typei,iscutmov,fraction_type},1)
                    a= Index_in_region{typei,sessioni,iscutmov,fraction_type}{regioni,:};
                    Response_amp_region_allsession(:,:,sessioni,regioni,typei)=squeeze(mean(A_r(:,:,a),3,'omitnan'));
                    Response_sd_region_allsession(:,:,sessioni,regioni,typei)=squeeze(std(A_r(:,:,a),[],3,'omitnan'));
                    Response_amp2_region_allsession(:,:,sessioni,regioni,typei)=squeeze(mean(mean(A_r(31:41,:,a),1,'omitnan'),3,'omitnan'));
                end
            end
            for align_typei=1
                for session_align=1:size(align_win,2)
                    switch align_typei
                        case 1
                            a=align_win(:,session_align);
                        case 2
                            a=align_win_go(:,session_align);
                        case 3
                            a=align_win_nogo(:,session_align);
                    end
                    a(isnan(a))=[];
                    Response_amp_align_allsession1(:,:,session_align,typei,align_typei)=mean(Response_amp_allsession(:,:,a,typei),3,'omitnan');
                    Response_sd_align_allsession1(:,:,session_align,typei,align_typei)=mean(Response_sd_allsession(:,:,a,typei),3,'omitnan');
                    Response_amp_region_allsession1(:,:,session_align,:,typei,align_typei)=mean(Response_amp_region_allsession(:,:,a,:,typei),3,'omitnan');
                    Response_sd_region_allsession1(:,:,session_align,:,typei,align_typei)=mean(Response_sd_region_allsession(:,:,a,:,typei),3,'omitnan');
                end
                for session_align=1:size(align_win,2)
                    switch align_typei
                        case 1
                            a=align_win(:,session_align);
                        case 2
                            a=align_win_go(:,session_align);
                        case 3
                            a=align_win_nogo(:,session_align);
                    end
                    a(isnan(a))=[];
                    Response_amp_align_allsession(:,session_align,:,typei,align_typei)=mean(Response_amp_align_allsession1(:,a,:,typei,align_typei),2,'omitnan');
                    Response_sd_align_allsession(:,session_align,:,typei,align_typei)=mean(Response_sd_align_allsession1(:,a,:,typei,align_typei),2,'omitnan');
                    Response_amp_region_align_allsession(:,session_align,:,:,typei,align_typei)=mean(Response_amp_region_allsession1(:,a,:,:,typei,align_typei),2,'omitnan');
                    Response_sd_region_align_allsession(:,session_align,:,:,typei,align_typei)=mean(Response_sd_region_allsession1(:,a,:,:,typei,align_typei),2,'omitnan');
                    
                end
            end
        end 
        a=Response_amp_all_session_align_all;a=cat(6,a,Response_amp_align_allsession);Response_amp_all_session_align_all=a;
        a=Response_sd_all_session_align_all;a=cat(6,a,Response_sd_align_allsession);Response_sd_all_session_align_all=a;
        a=Response_amp_region_all_session_align_all;a=cat(7,a,Response_amp_region_align_allsession);Response_amp_region_all_session_align_all=a;
        a=Response_sd_region_all_session_align_all;a=cat(7,a,Response_sd_region_align_allsession);Response_sd_region_all_session_align_all=a;
    end
end
% save([savepath,'/Summary_all_CRUR_singletrial.mat'],'A','path_all','Path','iscutmov','fraction_type','-v7.3');
save([savepath,'/Summary_all_CRUR_singletrial_append.mat'],'Response_amp_all_session_align_all','Response_sd_all_session_align_all','Response_amp_region_all_session_align_all','Response_sd_region_all_session_align_all');


%% plot typica trace
labels_typei={'CSon-avtivation','CSon-inhibition','CSoff-avtivation','CSoff-inhibition',...
    'US-activation','US-activation-cuttrail',...
    'CSon-up regulate','CSon-down regulate','CSon-stable regulate',...
    'CSoff-up regulate','CSoff-down regulate','CSoff-stable regulate',...
    'US-up regulate','US-down regulate','US-stable regulate',...,
    'US-up regulate-cuttrail','US-down regulate-cuttrail','US-stable regulate-cuttrail',...,
    'CSon-up regulate-all','CSon-down regulate-all','CSon-stable regulate-all',...
    'CSoff-up regulate-all','CSoff-down regulate-all','CSoff-stable regulate-all',...
    'US-up regulate-all','US-down regulate-all','US-stable regulate-all',...,
    'US-up regulate-cuttrail-all','US-down regulate-cuttrail-all','US-stable regulate-cuttrail-all',...
    'CSUS-shift','CSUS-shift-cuttrail','CSUS-shift-all','CSUS-shift-cuttrail-all'};
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
group_label={'Learner (n=8)','Faded-Learner (n=6)','Non-Learner (n=9)','Control (n=8)'};
align_label={'All trials','Acting trials','Non-acting trials'};
mean1= @(x)(mean(x,ndims(x),'omitnan'));
error1= @(x)(std(x,[],ndims(x),'omitnan')./sqrt(ndims(x)));
savepath_eps=checkpath('X:\calcium data 20230224\2');
savetype='jpg';
data=Response_amp_all_session_align_all;size(data)
clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);

for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 60];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
            name=[labels_typei{typei},'_',align_label{align_typei}];
            h=figure('position',[33,50,1800,1800],'name',name,'color','k');p=[];
            t=tiledlayout(length(ss),length(sessionx));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                for sessionj=1:length(sessionx)
                    dataL=squeeze(mean1(data(:,sessionj,ss(sessioni),typei,align_typei)));
                    sd=squeeze(error1(data(:,sessionj,ss(sessioni),typei,align_typei)));
                    x=[1:size(dataL,1)]*fs.ca;
                    a=[-1 2];
                    nexttile;patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(a) min(a) max(a)  max(a)],'b','Facealpha',0.4,'edgealpha',0);
                    shadedErrorBar(x,dataL,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                    p(sessioni)=plot(x,dataL,'color','r','linewidth',4);hold on;
                    title(sessionx(sessionj));
                    set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
                    ylim(a);xlim([15,60]*fs.ca);         
                end
            end
            %legend(p,sessionx(ss),'color','w','location','NorthEastOutside','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none');
            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(0.8,0.2,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name]);
            close(h)
    end
end

%%
if count_spatial_p
    path=Path{batchi};loc={};all_loc={};ID=[];
    for ii=1:length(path)
        nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)]
        id=find(strcmp( M.label,nn));
        loc{ii}=supervoxel(id,1:3);
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        all_loc{ii}=supervolxeli(:,1:3);
        ID(ii,1:length(id))=id;
    end
    [P_z,H_z,P_t,H_t,Dr_ij,Dn_ij]=statistic_P_in_spatical_loc(loc,all_loc);
    save([checkpath(fullfile(savepath,'p_spatial-types\')),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_spatial.mat'],'loc','all_loc','ID', 'Dn_ij','Dr_ij', 'P_z','H_z','P_t','H_t', '-v7.3');
end
%% demo
num=9;res=[0.66,0.66,10];v=[90,90];radium=floor(15/0.66);
load(fullfile(path_all{num},'para.mat'));load(fullfile(path_all{num},'env.mat'));
load(fullfile(path_all{num}, 'activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo')
load(fullfile(path_all{num},'singletrial','CR_ind_summary_singletrial.mat'),'CR_ind_up');
supervoxel_raw=env.supervoxel;
for ii=unique(supervoxel_raw(:,3))'
    id=find(supervoxel_raw(:,3)==ii);
    supervoxel_raw(id,3)=(ii-1)*10/0.66+1;
end
supervoxel_raw(:,1)=env.width-supervoxel_raw(:,1);
for ii=1:size(align_win,2)
    a=align_win(:,ii);a(isnan(a))=[];
    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((length(a)/3)),3);t.TileSpacing = 'compact';t.Padding = 'compact';
    for jj=1:length(a)
    ind=find(CR_ind_up(:,a(jj),iscutmov)==1);leg=[num2str(a(jj))];
    nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
    set(gca,'fontsize',16)
    end
end
%%
clc;
clear all;
iscutmov=1;fraction_type=1;
savepath='X:\calcium data 20230224\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load([savepath,'/Path.mat'])
load([savepath,'/Summary_all_CRUR_singletrial.mat'])
%% select fish
numL=[1,3,4,6,7,8,9,10];
numC=[12,13,14,16,17,18,19,20];numNL=[21:29];numFL=[30:35];
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
group_label={'L','FL','NL','C'};
align_label={'All','Go','NoGo'};
labels_typei_ADJ={'CS_avtivation','CS_inhibition','US_activation','CSupregulate','CSdownregulate','CSstableregulate','USupregulate','USdownregulate','USstableregulate','CSnewemergedneuron'};
load(fullfile(path{1},'CR_ind_summary_align.mat'),'Label_region');
regionxx=categorical(Label_region{1,1}(:,1));regionxx = reordercats(regionxx,Label_region{1,1}(:,1));
mean1= @(x)(mean(x,ndims(x),'omitnan'));
error1= @(x)(std(x,[],ndims(x),'omitnan')./sqrt(length(x)));
error2= @(x)(mean(x,ndims(x),'omitnan')-std(x,[],ndims(x),'omitnan'));
clr_group=[colorplus(389);colorplus(253);colorplus(233);colorplus(448)];    
clr_cmap =addcolorplus(301); addcolorplus(336);%cmap = GetColormap('jet',81);
clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
clr_CStypes = addcolorplus(275);%clr_CStypes= ColorMap(clr_CStypes,3);clr_CStypes = flipud(clr_CStypes);
clr_UStypes = addcolorplus(272);%clr_UStypes= ColorMap(clr_UStypes,3);clr_UStypes = flipud(clr_UStypes);
%%