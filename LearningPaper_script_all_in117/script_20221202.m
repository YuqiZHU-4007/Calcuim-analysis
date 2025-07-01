clc;clear all;close all;
res=[0.66,0.66,10];radium=floor(30);
load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
load('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\Path.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
issort=1;isplot=1;
for batchi=1;length(Path);
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi};nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
        if exist(fullfile(path_fishi,'brain_region_related_statistic.mat'),'file')
            %% load
            load(fullfile(path_fishi,'para.mat'));
            load(fullfile(path_fishi,'env.mat'));
            load(fullfile(path_fishi,'brain_region_related_statistic.mat'));
            supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
            %% para
            frame_ind_hab=frame.cs_start:frame.cs_end;frame_ind_test=frame.cs_start:frame.cs_end;frame_ind_acq=frame.cs_start:frame.us_start-1;
            stimCS=zeros(1,frame.per_cycle*trial.total(1));
            for ii=1:size(stimCS,2)/frame.per_cycle
                stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=0;%23
            end
            stimUS=ones(1,frame.per_cycle*trial.acq(1))*3;
            for ii=trial.acq(2):trial.acq(3)
                stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=0;%16
            end
            supervoxel_plot1=[env.supervoxel(:,2) env.height-env.supervoxel(:,1) (env.depth-env.supervoxel(:,3))];
            supervoxel_plot2=[(env.supervoxel(:,2))*res(2) (env.height-env.supervoxel(:,1))*res(1) (env.depth-env.supervoxel(:,3))*res(3)];
            supervoxel_plot11=[env.supervoxel(:,2) env.supervoxel(:,1) (env.supervoxel(:,3))];
            supervoxel_plot21=[(env.supervoxel(:,2))*res(2) (env.supervoxel(:,1))*res(1) (env.supervoxel(:,3))*res(3)];
            supervoxel_plot3=[floor(supervolxeli(:,2)/res(2)) floor(supervolxeli(:,1)/res(1)) floor(supervolxeli(:,3)/res(3)+1)];
            region_id=[5,6];
            for regioni=1:size(region_id,1)
                a=region_id(regioni,:); a(find(a==0))=[];
                if length(a)>1
                    Region=strcat(Label{region_id(regioni,1)},'-',Label{region_id(regioni,2)});
                else
                    Region=Label{region_id(regioni,1)};
                end
                disp(strcat('Runnning...',path_fishi,'...',Region));
                %!!!!!! savepath
                savepath=checkpath(fullfile(path_fishi,'clust_durCond',Region));
                load([savepath,'/clust_results.mat']);
                cluster_type={'PCA-hireatch','PCA-k-means','k-means','auto cluster'};
                cluster_typei=3;
                savepathh=fullfile(savepath,cluster_type{cluster_typei})
                load([savepathh,'/clust_results.mat']);
                clrmap = GetColormap(clrmap_name,max(gIX));
                %% sort gIX
                if issort
                    m=[];k=[];
                    for kk=unique(gIX)'
                        ind=find(gIX==kk);
                        m(:,kk)=mean(act_all_acq_block_CR_AUC(:,cIX(ind)),2);
                        pp = polyfit([1:size(act_all_acq_block_CR_AUC,1)]',m(:,kk),1);
                        k(:,kk)=pp(1);
                    end
                    [K,I]=sort(k);%clrmap_s=clrmap(I,:);
                    gIX_s=[];cIX_s=[];kk=1;
                    for iiii=I
                        ind=find(gIX==iiii);
                        gIX_s=[gIX_s;kk*ones(size(ind))];
                        cIX_s=[cIX_s;cIX(ind)'];
                        kk=kk+1;
                    end
                    %         [h,ratio]=plot_test_2(act_all_acq_block,cIX,gIX,frame,clrmap,4,[-0.02 0.02],true,true);
                    %         [h,ratio]=plot_test_2(act_all_acq_block,cIX_s,gIX_s,frame,clrmap,4,[-2 2],true,true);
                    %[h]=pushbutton_popupplot_Callback(M_cluster,cIX,gIX,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
                    %[h]=pushbutton_popupplot_Callback(M_cluster,cIX_s,gIX_s,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
                else
                    cIX_s=cIX;gIX_s=gIX;
                end
                %% plot
                cIX=cIX_s;gIX=gIX_s;savepath_figure=checkpath(fullfile(savepathh,'/'));
                %% trace of clusteri
                savepath_figure2=checkpath(fullfile(savepathh,'Heatmap & Trace dur Cond for each cluster'));
                savepath_figure3=checkpath(fullfile(savepathh,'Heatmap & Trace Pre and Post Cond for each cluster'));
                savepath_figure4=checkpath(fullfile(savepathh,'AUC for each cluster'));
                for jj=[12]
                    ind=find(gIX==jj);
                    figname=['Cluster ', num2str(jj)];
                    if exist([savepath_figure2, figname,'.png'],'file')==0
                        h=plot_avgtrace_heatmap(M_cluster(cIX(ind),:),act_all_acq_block(:,cIX(ind)),frame,frame_ind_acq,trial_ind_acq,startpoint,fs,frameb);
                        title(figname); %saveas(h,[savepath_figure2,figname],'png');%savefig(h,[savepath_figure2,figname]);
                        %close(h)
                    end
                    
                    
                    
                    
                    figname=['Cluster ', num2str(jj)];
                    if exist([savepath_figure3, figname,'.png'],'file')==0
                        h=plot_avgtrace_heatmap2(M_cluster_hab(cIX(ind),:),act_all_hab(:,cIX(ind)),frame,frame_ind_hab,trial_ind_hab,startpoint,fs,frameb,...
                            M_cluster_test(cIX(ind),:),act_all_test(:,cIX(ind)),frame_ind_test,trial_ind_test);
                        title(figname); %saveas(h,[savepath_figure3,figname],'png');%savefig(h,[savepath_figure3,figname]);
                        %close(h)
                        
                    end
                    figname=['Cluster ', num2str(jj)];
                    if exist([savepath_figure4, figname,'.png'],'file')==0
                        act_all_acq_block_gIX=act_all_acq_block(:,cIX(ind));
                        act_all_hab_gIX=act_all_hab(:,cIX(ind));
                        act_all_test_gIX=act_all_test(:,cIX(ind));
                        act_all_acq_block_gIX=reshape(act_all_acq_block_gIX,frame.per_cycle,trial.acq(1),[]);
                        act_all_hab_gIX=reshape(act_all_hab_gIX,frame.per_cycle,trial.hab(1),[]);
                        act_all_test_gIX=reshape(act_all_test_gIX,frame.per_cycle,trial.test(1),[]);
                        act_all=cat(2,act_all_hab_gIX,act_all_acq_block_gIX,act_all_test_gIX);
                        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
                        [h,AUC,AUC_go,AUC_nogo,p_AUC,p_AUC_go,p_AUC_nogo]=plot_AUC(act_all,CR,trial,frame);
                        saveas(h,[savepath_figure4,figname],'png');
                        save([savepath_figure4,figname,'.mat'],'AUC','AUC_go','AUC_nogo','p_AUC','p_AUC_go','p_AUC_nogo', '-v7.3');
                    end
                end
            end
        end
    end
end
