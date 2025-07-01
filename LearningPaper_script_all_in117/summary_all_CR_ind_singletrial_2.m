
clc;
clear all;
%% load all 
iscutmov=1;fraction_type=1;
savepath='X:\calcium data 20230224\';
%savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
 warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
 %warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';

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
Response_amp_all=[];Response_amp_region_all=[];Response_amp_align_all=[];Response_amp_region_align_all=[];
Response_amp_region_smooth_all=nan(size(Response_amp_region_all));Response_amp_region_smooth_align_all=[];
start_index_all=[];start_index_align_all=[];
% A=setfield(A,'NUM_CR_UR',[]);A=setfield(A,'fraction_CR_UR',[]);
% A=setfield(A,'NUM_CR_UR_align',[]);A=setfield(A,'fraction_CR_UR_align',[]);
% A=setfield(A,'CS_response_amp',[]);A=setfield(A,'CS_response_amp_align',[]);
% A=setfield(A,'region_CS_response_amp',[]);A=setfield(A,'region_CS_response_amp_align',[]);
% A=setfield(A,'supervol_type',[]);A=setfield(A,'supervol_type_align',[]);A=setfield(A,'region_supervol',[]);A=setfield(A,'region_supervol_align',[]);
% A=setfield(A,'supervoltemp_type',[]);A=setfield(A,'supervoltemp_type_align',[]);A=setfield(A,'region_supervoltemp',[]);A=setfield(A,'region_supervoltemp_align',[]);
for groupi=[1,2,3,4]
    path=Path{groupi};
    for fishi=1:length(path)
        p=path{fishi};path_all{kk,1}=p;kk=kk+1;
        load(fullfile(p,'singletrial','CR_ind_summary_singletrial.mat'),'A_r','ind_all_emergeCSUS_inthissession','Index_in_region_emerged_inthissession','Fraction_in_region_type_emerged_inthissession','US_activation_response_cuttrail', 'CSon_activation_response');
        load(fullfile(p, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo');
        load(fullfile(p, '/env.mat'));supervol=env.supervoxel(:,1:3);
        nn=[p(end-14:end-7),p(end-5:end-1)];
        supervol_temp=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        NUM_CR_UR=[];fraction_CR_UR=[];NUM_CR_UR_align=[];fraction_CR_UR_align=[];
        Response_amp=[];Response_amp_region=nan(frame.per_cycle,trial.total,30,size(ind_all_emergeCSUS_inthissession,1));
        Response_amp_align=[];Response_amp_region_align=[];
        CS_response_amp=[];CS_response_amp_align=[];region_CS_response_amp=[];region_CS_response_amp_align=[];
        supervol_type={};supervoltemp_type={};region_supervol={};region_supervoltemp={};
        supervol_type_align={};supervoltemp_type_align={};region_supervol_align={};region_supervoltemp_align={};
        start_index=nan(3,30,trial.total,size(ind_all_emergeCSUS_inthissession,1));Response_amp_region_smooth=[];start_index_align=[];Response_amp_region_smooth_align=[];
        for typei=1:size(ind_all_emergeCSUS_inthissession,1)
            for sessioni=1:size(Fraction_in_region_type_emerged_inthissession{typei,iscutmov,fraction_type,1},2)
                NUM_CR_UR(sessioni,typei)=length(ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1})./size(supervol,1);
                fraction_CR_UR(sessioni,:,typei)=Fraction_in_region_type_emerged_inthissession{typei,iscutmov,fraction_type,1}(:,sessioni);
                CS_response_amp(sessioni,typei)=squeeze(mean( mean(CSon_activation_response{sessioni,iscutmov}(:,:,ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1}),1,'omitnan'),3,'omitnan'));
                supervol_type{typei,sessioni}=supervol(ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1},:);
                supervoltemp_type{typei,sessioni}=supervol_temp(ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1},:);
                Response_amp(:,sessioni,typei)=squeeze(mean(A_r(:,sessioni,ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1}),3,'omitnan'));
                for regioni=1:size(Fraction_in_region_type_emerged_inthissession{typei,iscutmov,fraction_type,1},1)
                    a= Index_in_region_emerged_inthissession{typei,sessioni,iscutmov,fraction_type,1}{regioni,:};
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
                    for regioni=1:size(region_supervol,2)
                        region_supervol_align{typei,regioni,session_align,align_typei}={region_supervol{typei,regioni,a}};
                        region_supervoltemp_align{typei,regioni,session_align,align_typei}={region_supervoltemp{typei,regioni,a}};
                        Response_amp_region_align(:,session_align,regioni,typei,align_typei)=mean(Response_amp_region(:,a,regioni,typei),2,'omitnan');
                        Response_amp_region_smooth_align(:,session_align,regioni,typei,align_typei)=mean(Response_amp_region_smooth(:,a,regioni,typei),2,'omitnan');
                    end
                    start_index_align(:,:,session_align,typei,align_typei)=mean(start_index(:,:,a,typei),3,'omitnan');
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
save([savepath,'/Summary_all_CRUR_singletrial_ind_all_emergeCSUS_inthissession.mat'],'Path','path_all','iscutmov','fraction_type','base_win','onset_win',...
    'NUM_CR_UR_all','NUM_CR_UR_align_all','fraction_CR_UR_all','fraction_CR_UR_align_all',...
'CS_response_amp_all','CS_response_amp_align_all','region_CS_response_amp_all','region_CS_response_amp_align_all',...
'region_supervol_all','region_supervol_align_all',...
'Response_amp_all','Response_amp_region_all','Response_amp_align_all','Response_amp_region_align_all','-v7.3');
save([savepath,'/Summary_all_CRUR_singletrial_ind_all_emergeCSUS_inthissession.mat'],'Path','path_all','iscutmov','fraction_type','base_win','onset_win',...
    'NUM_CR_UR_all','NUM_CR_UR_align_all','fraction_CR_UR_all','fraction_CR_UR_align_all',...
'CS_response_amp_all','CS_response_amp_align_all','region_CS_response_amp_all','region_CS_response_amp_align_all',...
'supervol_type_all','supervol_type_align_all','region_supervol_all','region_supervol_align_all',...
'supervoltemp_type_all','supervoltemp_type_align_all','region_supervoltemp_all','region_supervoltemp_align_all',...
'Response_amp_all','Response_amp_region_all','Response_amp_align_all','Response_amp_region_align_all',...
'Response_amp_region_smooth_all','Response_amp_region_smooth_align_all','start_index_all','start_index_align_all','-v7.3');


clc;
clear all;
%% load all
iscutmov=1;fraction_type=1;
%savepath='X:\calcium data 20230224\';%
savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
%warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load(fullfile(Path{1}{1},'para.mat'));
Response_amp_all_session_align_all=[];

for groupi=1:4
    path=Path{groupi};
    for fishi=1:length(path)
        p=path{fishi};
        load(fullfile(p,'singletrial','CR_ind_summary_singletrial.mat'),'A_r','ind_all_emergeCSUS_inthissession','Index_in_region_emerged_inthissession','Fraction_in_region_type_emerged_inthissession');
        load(fullfile(p, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo');
        load(fullfile(p, '/env.mat'));supervol=env.supervoxel(:,1:3);
        nn=[p(end-14:end-7),p(end-5:end-1)];
        Response_amp_allsession=[];Response_amp_align_allsession=[]; Response_amp_align_allsession1=[];Response_sd_align_allsession1=[];
        Response_sd_allsession=[];Response_sd_align_allsession=[];
        for typei=7;%1:size(ind_all_emergeCSUS_inthissession,1)
            for sessioni=1:size(Fraction_in_region_type_emerged_inthissession{typei,iscutmov,fraction_type},2)
                Response_amp_allsession(:,:,sessioni,typei)=squeeze(mean(A_r(:,:,ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1}),3,'omitnan'));     
                Response_sd_allsession(:,:,sessioni,typei)=squeeze(std(A_r(:,:,ind_all_emergeCSUS_inthissession{typei,sessioni,iscutmov,1}),[],3,'omitnan'));     

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
                    Response_amp_align_allsession1(:,:,session_align,typei,align_typei)=mean(Response_amp_allsession(:,:,a,typei),3,'omitnan');
                    Response_sd_align_allsession1(:,:,session_align,typei,align_typei)=mean(Response_sd_allsession(:,:,a,typei),3,'omitnan');

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
                end
            end
        end 
        a=Response_amp_all_session_align_all;a=cat(6,a,Response_amp_align_allsession);Response_amp_all_session_align_all=a;
    end
end
% save([savepath,'/Summary_all_CRUR_singletrial.mat'],'A','path_all','Path','iscutmov','fraction_type','-v7.3');
save([savepath,'/Summary_all_CRUR_singletrial_ind_all_emergeCSUS_inthissession_append.mat'],'Response_amp_all_session_align_all','-v7.3');
