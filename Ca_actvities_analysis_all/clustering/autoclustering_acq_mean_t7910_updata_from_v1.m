%AutoClustering of all learner's CR during conditioning
%refer to Chen,2018,Neuron
%auto clustering in T7910
%other related file was 'Get_all_fish_act.mat' & 'Get_all_fish_act_control.m'
addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
clc;close all;clear all;
set(0,'defaultfigurecolor','w');
issort=true;ismerge=false;merge_matrix=[];
actbatch={};envbatch={};behavbatch={};parabatch={};
% actbatch=batch_code('txt of <activties_new> path') ;
% envbatch=batch_code('txt of <location> path') ;
% behavbatch=batch_code('txt of <behavior> path') ;
% parabatch=batch_code('txt of <para> path') ;
load('K:\2.Fear conditioning_huc_dimming\autoclustering_results\path\path.mat');%load path.mat
savepath='K:\2.Fear conditioning_huc_dimming\autoclustering_results\';
%% correct actpath
name='activities_aft_process';
for ii=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{ii});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{ii}=fullfile(fpath,[fname,ext]);
    end
end
%% extra all acti
colnum=3000;
cut_move=true;
is_eatra_act_all=false;
act_all=[];index_all=[];
act_all_hab=[];act_all_hab_mean=[];
act_all_acq=[];act_all_acq_block_mean=[];
act_all_tst=[];act_all_tst_mean=[];
env_all.supervoxel=[];
area_all=struct;
for kk=2:size(actbatch,1)
    load(actbatch{kk});
    load(behavbatch{kk});
    load(parabatch{kk});
    act=activities_preCS_dfdf_aftcorrect';
    index_all=[index_all, ones(1,size(act',2))*kk];
    load(envbatch{kk});
    env_all.supervoxel=[env_all.supervoxel; env.supervoxel];
    if is_eatra_act_all
        if size(act,2)>colnum %size(act,2)~=size(act_all,1) && ~isempty(act_all)
            act(:,colnum+1:end)=[];
            warning([actbatch{kk} ' acti was cut into ' num2str(size(act_all,1))]);
        end
        act_all=[act_all , act'];
    end
    act=activities_preCS_dfdf_aftcorrect';
    i_ind=3;
    ind_cut_trial_preCS=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
    ind_cut_trial_preCS=unique(ind_cut_trial_preCS);
    ind_cut_trial_CS=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start) & re_startpoint(:,2)<(frameb.us_start-1)),1);%%%%%%%%%%%%%%%%
    ind_cut_trial_CS=unique(ind_cut_trial_CS);
    ind_cut_trial=union(ind_cut_trial_preCS,ind_cut_trial_CS);
    [rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
    %     if ~exist('rawacti') && i_ind>1 && ~exist('area')
    %         [rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial_preCS);
    %     end
    for ii=1:i_ind
        switch ii
            case 1
                trial_ind=trial.hab(2):min(trial.hab(2)+6-1,trial.hab(3)); %%trial.test(2):trial.test(2)+2
                act_CS_hab= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,[],size(act,1));
                %act_all_hab(:,:,end+1:end+size(act_CS_hab,3))=act_CS_hab;
                act_all_hab=cat(3,act_all_hab,act_CS_hab);
                act_CS_hab_mean=reshape(mean(act_CS_hab,2),frame.per_cycle,[]);
                act_all_hab_mean=[act_all_hab_mean,act_CS_hab_mean];
            case 2
                trial_ind=trial.acq(2):trial.acq(3); %%trial.test(2):trial.test(2)+2
                aa=[];
                for jj=trial_ind
                    %a=act(:,(jj-1)*frame.per_cycle+frame.cs_start:(jj-1)*frame.per_cycle+frame.us_start-1)';
                    a=act(:,(jj-1)*frame.per_cycle+1:(jj)*frame.per_cycle)';
                    aa=[aa;a];
                end
                act_all_acq=[act_all_acq aa];
                for jj=1:size(act',2)
                    act_all_acq_block_mean(:,:,size(act_all_acq_block_mean,3)+1)=rawacti.acq_mean_blocks{jj};
                end
                name=fieldnames(area);
                for jj=1:length(name)
                    if ~isfield(area_all,name{jj})
                        area_all=setfield(area_all,name{jj},[]);
                    end
                    a=getfield(area,name{jj});
                    aa=getfield(area_all,name{jj});
                    area_all=setfield(area_all,name{jj},[aa a]);
                end
            case 3
                trial_ind=trial.test(2):min(trial.test(2)+3-1,trial.test(3)); %%trial.test(2):trial.test(2)+2
                act_CS_hab= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,[],size(act,1));
                act_all_tst(:,:,end+1:end+size(act_CS_hab,3))=act_CS_hab;
                act_CS_hab_mean=reshape(mean(act_CS_hab,2),frame.per_cycle,[]);
                act_all_tst_mean=[act_all_tst_mean,act_CS_hab_mean];
        end
        %act_CS_hab_mean=zeros(frame.per_cycle,size(act,1));
    end
end
act_all_acq_block_mean(:,:,1)=[];act_all_tst(:,:,1)=[];
save([savepath '\act'],'act_all','index_all','act_all_hab','act_all_hab_mean','act_all_acq','act_all_acq_block_mean','act_all_tst','act_all_tst_mean','env_all','area_all','-v7.3');
%save(['F:\DUlab\Progress report\Progress Report201905-201906\' '\act'],'act_all','index','act_all_hab','act_all_hab_mean','act_all_acq','act_all_acq_block_mean','act_all_tst','act_all_tst_mean','env_all','area_all','-v7.3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classification
ref_win=ceil(0/fs.ca+1):frame.cs_start-1;%取时间段的baseline作为参考判断event
time=([1:frame.per_cycle]-frame.cs_start)*fs.ca;
area_win_hab=frame.cs_start:frame.cs_end-1;
stimCS=zeros(1,frame.per_cycle*trial.total(1));
for ii=1:size(stimCS,2)/frame.per_cycle
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=0;%23
end
stimUS=ones(1,frame.per_cycle*trial.acq(1))*3;
for ii=trial.acq(2):trial.acq(3)
    stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=0;%16
end
%%%%%Acq.
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
%act_all_acq_block_mean(:,:,1)=[];
M_0=reshape(act_all_acq_block_mean(frame_ind,:,:),length(frame_ind)*size(act_all_acq_block_mean,2),[])';
M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
%%%clustering
numK=[];
ind=index_all;find(index_all==2);
ind=ind([1:50:length(ind)]);
[gIX,~,numK,~]=find_best_k_in_range(M_0(ind,:),1:40);
cIX=ind;cIX_reg = (1:size(M_0,1))';%ind;
if isempty(numK)
    numK=20;
end
masterthres=0.6;
para=struct;
para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
%cIX=(1:size(M_0,1))';gIX=(1:size(M_0,1))'; cIX_reg = (1:size(M_0,1))';
for ii=1:5
    %[cIX,gIX] = AutoClustering(randi(size(M_0,1),[1,50000]),[1:50000],M_0,cIX_reg,1,para,1,0.7);
    [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
    cIX_reg = (1:size(M_0,1))';
end
result.cIX=cIX;result.gIX=gIX;%!!!!!!!!!!!!!!!!!!!!!!!!

%% sort gIX
if issort
    A=area_all.CS_acq_block;m=[];
    for kk=unique(result.gIX)'
        ind=find(result.gIX==kk);
        m(:,kk)=mean(A(:,result.cIX(ind)),2);
        p = polyfit([1:size(A,1)]',m(:,kk),1);
        k(:,kk)=p(1);
    end
    
    [K,I]=sort(k);%clrmap_s=clrmap(I,:);
    gIX_s=[];cIX_s=[];kk=1;
    for ii=I
        ind=find(result.gIX==ii);
        gIX_s=[gIX_s;kk*ones(size(ind))];
        cIX_s=[cIX_s;result.cIX(ind)];
        kk=kk+1;
    end
    gIX=gIX_s;cIX=cIX_s;
    % figure,a=[-0.01 0.1];plot_test_1([1:5],A,cIX,gIX,clrmap,a,frame,fs.ca);
    % figure,a=[-0.01 0.1];plot_test_1([1:5],A,cIX_s,gIX_s,clrmap_s,a,frame,fs.ca);
    % [h]=pushbutton_popupplot_Callback(B,cIX_s,gIX_s,clrmap_s,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
end
%% merge gIX
if ismerge
    merge_matrix
    gIX_m=gIX;
    for ii=1:size(merge_matrix,1)
        merge_ind=merge_matrix(ii,:);
        for jj=1:length(merge_ind)
            ind=find(gIX==merge_ind(jj));
            if ~isempty(ind)
                gIX_m(ind)=ii;
            end
        end
    end
    gIX=gIX_m;
end
%% some plot
%cIX=result.cIX;gIX=result.gIX;%!!!!!!!!!!!!!!!!!!!!!!!!
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX));
ind_control=setdiff(1:length(index_all),cIX)';%length(ind)+length(cIX)== length(index_all)
clr_control=[0,0,0];
index_all=index_all';
cIX_all=[cIX;ind_control];
gIX_all=[gIX;(max(gIX)+1)*ones(size(ind_control))];
clrmap_all=cat(1,clrmap,clr_control);
length(unique(gIX))
%cIX(35670:35679)=[];gIX(35670:35679)=[];
vi='off';
for iii=1:2
    switch iii
        case 1
            cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
        case 2
            cIX_iii=cIX_all;gIX_iii=gIX_all;clrmap_iii=clrmap_all;
        case 3
            cIX_iii=cIX;gIX_iii=gIX;
            cIX_iii(find(gIX==11 | gIX==12 | gIX==13))=[];gIX_iii(find(gIX==11 | gIX==12 | gIX==13))=[];clrmap_iii =clrmap;
            unique(gIX_iii)
    end
    outputpath=checkpath([savepath '\figures_masterthres' num2str(masterthres) '_' num2str(iii,'%02d')]);
    %1
    [h,~, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M_0(cIX_iii,:));saveas(h,[outputpath '\f1'],'fig');
    %2
    B=act_all_acq';
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[outputpath '\f2'],'fig');
    %3
    h=figure;B=cat(1,area_all.CS_hab_tst(1,:),area_all.CS_acq_block,area_all.CS_hab_tst(2,:));%B=[area.CS_hab_tst(1,:); area.CS_acq_block ;area.CS_hab_tst(2,:)];
    a=[-0.01 0.1];plot_test_1([1:size(B,1)],B,cIX_iii,gIX_iii,clrmap_iii,a,frame,fs.ca,2);saveas(h,[outputpath '\f3'],'fig');
    %4
    [h,ratio]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,1,[-0.02 0.04],1,1);saveas(h,[outputpath '\f41'],'fig');
    [h,ratio]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.02 0.04],1,1);saveas(h,[outputpath '\f42'],'fig');
    %5
    [h,~]=plot_test_2(act_all_hab_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.02 0.04],1,1);saveas(h,[outputpath '\f5'],'fig');
    %6
    [h,~]=plot_test_2(act_all_tst_mean,cIX_iii,gIX_iii,frame,clrmap_iii,2,[-0.02 0.04],1,1);saveas(h,[outputpath '\f6'],'fig');
    %7
    [h1,h2,ratio_fish]=plot_test_3(act_all_acq,cIX_iii,gIX_iii,clrmap_iii,index_all,envbatch,env_all,outputpath,fs,stimCS,stimUS);
end
save([savepath '\clust_masterthres_' num2str(masterthres) '.mat'],'result','cIX','gIX','ratio','ratio_fish','para','masterthres');
%%%%%%%chose signle cluster num.
num_clust=[4,6,10,16];cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;a=[];
for kk=3%unique(index_all)
    [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
    load(envbatch{kk});
    act_trace=act_all_acq';
    ind_cls_in_this_fish=i_cIX;ind_add=[];
    for ii=num_clust
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
        %ind_cls_in_this_fish( find(gIX_iii(i_cIX)~=ii))=[];
    end
    ind_cls_in_this_fish=i_cIX(ind_add);
    colorind=unique(gIX_iii(ind_cls_in_this_fish));
    %a=a+length(ind_cls_in_this_fish);
    %a=[a;cIX_iii(ind_cls_in_this_fish)];
    %h=pushbutton_popupplot_Callback(act_trace,cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),clrmap_iii(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
    h=DrawTiledPics_zyq_20190530(cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),[1:size(act_trace,1)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap_iii);
    saveas(h,[savepath '\fish' num2str(kk,'%02d') '_' num2str(num_clust,'%02d')],'fig');
end

num_clust=[1:17];cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;a=[];
for ii=num_clust
    for kk=3%unique(index_all)
        [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
        load(envbatch{kk});
        act_trace=act_all_acq';
        ind_cls_in_this_fish=i_cIX;ind_add=[];
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
        %ind_cls_in_this_fish( find(gIX_iii(i_cIX)~=ii))=[];
        ind_cls_in_this_fish=i_cIX(ind_add);
        colorind=unique(gIX_iii(ind_cls_in_this_fish));
        %a=a+length(ind_cls_in_this_fish);
        %a=[a;cIX_iii(ind_cls_in_this_fish)];
        %h=pushbutton_popupplot_Callback(act_trace,cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),clrmap_iii(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
        h=DrawTiledPics_zyq_20190530(cIX_iii(ind_cls_in_this_fish),gIX_iii(ind_cls_in_this_fish),[1:size(act_trace,1)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap_iii);
        saveas(h,[savepath '\fish' num2str(kk,'%02d') '_' num2str(ii,'%02d')],'fig');
    end
end
