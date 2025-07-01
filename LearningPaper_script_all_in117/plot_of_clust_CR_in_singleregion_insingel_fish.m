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
            region_id=[3,4;5,6;25,0;26,0;7,8;13,0;29,0]%[3,4;5,6;25,0;26,0;7,8;13,0];%[3,4;5,6;7,8;9,10;11,12;13,0;14,0;15,0;16,0;17,0;21,0;22,0;23,24;25,0;26,0;29,0];%[3,4;5,6;25,0;26,0;7,8;13,0];%[29,0];%30,0;9,10;11,12;23,24;14,0;15,0;16,0;17,0;27,28;22,0;21,0];%
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
                if exist([savepath,'/clust_results.mat'],'file')~=0
                load([savepath,'/clust_results.mat']);
                cluster_type={'PCA-hireatch','PCA-k-means','k-means','auto cluster'};
                for cluster_typei=3;length(cluster_type);
                    savepathh=fullfile(savepath,cluster_type{cluster_typei})
                    if exist([savepathh,'/clust_results.mat'],'file')
                        load([savepathh,'/clust_results.mat']);
                        clrmap = GetColormap(clrmap_name,max(gIX));
                        if isplot==1
                            for isploti=1
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
%                                             [h,ratio]=plot_test_2(act_all_acq_block,cIX,gIX,frame,clrmap,4,[-0.02 0.02],true,true);
%                                             [h,ratio]=plot_test_2(act_all_acq_block,cIX_s,gIX_s,frame,clrmap,4,[-2 2],true,true);
%                                     [h]=pushbutton_popupplot_Callback(M_cluster,cIX,gIX,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
%                                     [h]=pushbutton_popupplot_Callback(M_cluster,cIX_s,gIX_s,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
                                else
                                    cIX_s=cIX;gIX_s=gIX;
                                end
                                save([savepathh,'/clust_results_s.mat'],'cIX_s','gIX_s');
                                
                                %% plot
                                cIX=cIX_s;gIX=gIX_s;savepath_figure=checkpath(fullfile(savepathh,'/'));
                                %% heatmap % mean trace
                                if exist([savepath_figure,'Heatmap & mean trace','.png'],'file') ==0              
                                h=plot_avgtrace_heatmap(M_cluster,act_all_acq_block,frame,frame_ind_acq,trial_ind_acq,startpoint,fs,frameb);
                                saveas(h,[savepath_figure,'Heatmap & mean trace'],'png');%savefig(h,[savepath_figure,'Heatmap & mean trace']); 
                                close(h)
                                end
                                %% t-sne结果:scatter of M_cluster_dr
                                if exist('M_cluster_dr','var')
                                    if exist([savepath_figure,'T-sne','.png'],'file')==0
                                    h=figure, gscatter(M_cluster_dr(:,1),M_cluster_dr(:,2),gIX);
                                    xlabel('Tsne-1');ylabel('Tsne-2');set(gca,'fontsize',16);legend off;
                                    saveas(h,[savepath_figure,'T-sne'],'png');%savefig(h,[savepath_figure,'T-sne']); 
                                    close(h)
                                    end
                                end
                                %% corrmap
                                figname='Corrmap';
                                if exist([savepath_figure, figname,'.png'],'file')==0
                                    h=figure('position',[484,403,400,350]);%subplot(1,2,1);imagesc(Corr_map(cIX,cIX),[min(Corr_map(:)) max(Corr_map(:))]);colormap('jet');colorbar;
                                    [idxAs,idxAa_ind]=sort(gIX);
                                    imagesc(Corr_map(cIX(idxAa_ind),cIX(idxAa_ind)),[min(Corr_map(:)) max(Corr_map(:))]);colormap('jet');colorbar;title(nn);
                                    saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]);
                                    close(h)
                                end
                                
                                %% imagesc
                                % y=roundn([min(M_cluster(:)) max(M_cluster(:))],1);
                                % h=figure,subplot(2,1,1),imagesc(M_cluster,y);colormap('hot');colorbar;
                                % [~,a]=sort(gIX);b=gIX(a);
                                % subplot(2,1,2),imagesc(M_cluster(cIX(a),:),y);colormap('hot');colorbar;
                                %% 占比
                                figname='Fraction of clustertypes';
                                if exist([savepath_figure, figname,'.png'],'file')==0
                                    h=figure;hh=histogram(gIX,'Normalization','probability','BinWidth',1,'BinEdge',[0.5:1:max(gIX)+0.5]);
                                    fraction=histcounts(gIX,'Normalization','probability','BinWidth',1,'BinEdge',[0.5:1:max(gIX)+0.5])*length(gIX);
                                    text(hh.BinEdges(1:end-1),hh.Values,num2str(fraction'),'HorizontalAlignment','center','VerticalAlignment','bottom');
                                    ylabel('Neurons(%)');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2,'xtick',1:max(gIX));box off;
                                    saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]);
                                    close(h)
                                    save([savepath_figure,figname,'.mat'],'fraction', '-v7.3');
                                end

                                %% mapback
                                figname='Mapback of all clusters to temp';
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                h=DrawTiledPics_zyq_20190530(neuron_id_in_regioni(cIX),gIX,[1:size(supervoxel_plot3,1)],supervoxel_plot3,temp_env.env.vol,clrmap);text(10,10,nn,'color','w');
                                saveas(h,[savepath_figure,figname],'png');savefig(h,[savepath_figure,figname]);
                                close(h)
                                 end
                                 figname='Mapback of all clusters to raw vol';   
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                     h=DrawTiledPics_zyq_20190530(neuron_id_in_regioni(cIX),gIX,[1:size(supervoxel_plot11,1)],supervoxel_plot11,env.vol,clrmap);text(10,10,nn,'color','w');
                                     saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]);
                                     close(h)
                                 end
                                
                                figname='Scatter Mapback of all clusters to raw vol';
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                h=plot_scatter_mapback_all([],supervoxel_plot2,neuron_id_in_regioni,cIX,gIX,clrmap,nn);
                                    saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]); 
                                    close(h)
                                 end
                                figname='Scatter Mapback of all clusters to temp';   
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                h=plot_scatter_mapback_all([],supervolxeli,neuron_id_in_regioni,cIX,gIX,clrmap,nn);
                                 saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]); 
                                 close(h)
                                 end
                                %% trace
                                figname='Trace of all clusters'; 
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                [h,ratio]=plot_test_2(act_all_acq_block,cIX,gIX,frame,clrmap,4,[-2 2],true,true);
                                   saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]); 
                                   close(h)
                                 end
                                
                                figname='Heatmap of all clusters'; 
                                 if exist([savepath_figure, figname,'.png'],'file')==0
                                [h]=pushbutton_popupplot_Callback(M_cluster,cIX,gIX,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);text(4028,7268,nn);
                                   saveas(h,[savepath_figure,figname],'png');%savefig(h,[savepath_figure,figname]); 
                                   close(h)
                                 end
                                
                                %% mapback to raw vol
                                savepath_figure2=checkpath(fullfile(savepathh,'Mapback to raw vol for each cluster','/'));
                                savepath_figure3=checkpath(fullfile(savepathh,'Scatter mapback to raw vol for each cluster','/'));
                                for jj=unique(gIX)'
                                    ind=find(gIX==jj);
                                    figname=['Cluster ', num2str(jj)];
                                    if exist([savepath_figure2, figname,'.png'],'file')==0
                                        h=DrawTiledPics_zyq_20190530(neuron_id_in_regioni(cIX(ind)),gIX(ind),[1:size(supervoxel_plot11,1)],supervoxel_plot11,env.vol,clrmap);text(10,10,nn,'color','w');
                                        saveas(h,[savepath_figure2,figname],'png');%savefig(h,[savepath_figure2,figname]);
                                        close(h)
                                    end
                                    figname=['Cluster ', num2str(jj)];
                                    if exist([savepath_figure3, figname,'.png'],'file')==0
                                        h=plot_scatter_mapback([],supervoxel_plot2,neuron_id_in_regioni,cIX,ind,['Cluster ',num2str(jj),' (n=',num2str(length(ind)),')']);
                                        saveas(h,[savepath_figure3,figname],'png');%savefig(h,[savepath_figure3,figname]);
                                        close(h)
                                    end
                                end
                                %% mapback to temp
                                X_Y_Z_Density={};
                                savepath_figure2=checkpath(fullfile(savepathh,'Mapback to temp for each cluster','/'));
                                savepath_figure3=checkpath(fullfile(savepathh,'Scatter mapback to temp for each cluster','/'));
                                savepath_figure4=checkpath(fullfile(savepathh,'Scatter mapback density to temp for each cluster','/'));
                                for jj=unique(gIX)'
                                    ind=find(gIX==jj);
                                    figname=['Cluster ', num2str(jj)];  
                                     if exist([savepath_figure2, figname,'.png'],'file')==0
                                    h=DrawTiledPics_zyq_20190530(neuron_id_in_regioni(cIX(ind)),gIX(ind),[1:size(supervoxel_plot3,1)],supervoxel_plot3,temp_env.env.vol,clrmap);text(10,10,nn,'color','w');
                                     saveas(h,[savepath_figure2,figname],'png');%savefig(h,[savepath_figure2,figname]); 
                                     close(h)
                                     end
                                    
                                     figname=['Cluster ', num2str(jj)]; 
                                      if exist([savepath_figure3, figname,'.png'],'file')==0
                                     h=plot_scatter_mapback([],supervolxeli,neuron_id_in_regioni,cIX,ind,['Cluster ',num2str(jj),' (n=',num2str(length(ind)),')']);
                                      saveas(h,[savepath_figure3,figname],'png');%savefig(h,[savepath_figure3,figname]); 
                                      close(h)
                                      end
                                    
                                      figname=['Cluster ', num2str(jj)];
                                      if exist([savepath_figure4, figname,'.png'],'file')==0
                                          nbin=size(temp_env.env.vol);h=figure;
                                          for type=1:3
                                              a=supervolxeli(neuron_id_in_regioni(cIX(ind)),type);
                                              [N,edges]=histcounts(a,ceil(nbin(type)/10),'BinWidth',10);
                                              x=edges(1:end-1)'; y=N';z=ones(length(x),1)*length(ind);
                                              X_Y_Z_Density{jj,type}=[x,y,z];
                                              subplot(3,1,type);plot(x,y,'Linewidth',1.5,'color','r');xlim([1,nbin(type)*res(type)])
                                          end
                                          saveas(h,[savepath_figure4,figname],'png');%savefig(h,[savepath_figure4,figname]);
                                          close(h)
                                          save([savepath_figure4,'density.mat'],'X_Y_Z_Density', '-v7.3');
                                      end
                                end
                                %% trace of clusteri
                                savepath_figure2=checkpath(fullfile(savepathh,'Heatmap & Trace dur Cond for each cluster'));
                                savepath_figure3=checkpath(fullfile(savepathh,'Heatmap & Trace Pre and Post Cond for each cluster'));
                                savepath_figure4=checkpath(fullfile(savepathh,'AUC for each cluster'));
                                for jj=unique(gIX)'
                                    ind=find(gIX==jj);
                                    figname=['Cluster ', num2str(jj)];
                                    if exist([savepath_figure2, figname,'.png'],'file')==0
                                        h=plot_avgtrace_heatmap(M_cluster(cIX(ind),:),act_all_acq_block(:,cIX(ind)),frame,frame_ind_acq,trial_ind_acq,startpoint,fs,frameb);
                                       title(figname); saveas(h,[savepath_figure2,figname],'png');%savefig(h,[savepath_figure2,figname]);
                                        close(h)
                                    end
                                    figname=['Cluster ', num2str(jj)];
                                    if exist([savepath_figure3, figname,'.png'],'file')==0
                                        h=plot_avgtrace_heatmap2(M_cluster_hab(cIX(ind),:),act_all_hab(:,cIX(ind)),frame,frame_ind_hab,trial_ind_hab,startpoint,fs,frameb,...
                                            M_cluster_test(cIX(ind),:),act_all_test(:,cIX(ind)),frame_ind_test,trial_ind_test);
                                        title(figname); saveas(h,[savepath_figure3,figname],'png');%savefig(h,[savepath_figure3,figname]);
                                        close(h)
                                    end
                                    figname=['Cluster ', num2str(jj)];
                                    if length(ind)>=2 %exist([savepath_figure4, figname,'.png'],'file')==0 && 
                                        act_all_acq_block_gIX=act_all_acq_block(:,cIX(ind));
                                        act_all_hab_gIX=act_all_hab(:,cIX(ind));
                                        act_all_test_gIX=act_all_test(:,cIX(ind));
                                        act_all_acq_block_gIX=reshape(act_all_acq_block_gIX,frame.per_cycle,trial.acq(1),[]);
                                        act_all_hab_gIX=reshape(act_all_hab_gIX,frame.per_cycle,trial.hab(1),[]);
                                        act_all_test_gIX=reshape(act_all_test_gIX,frame.per_cycle,trial.test(1),[]);
                                        act_all=cat(2,act_all_hab_gIX,act_all_acq_block_gIX,act_all_test_gIX);
                                        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
                                        [h,AUC,AUC_go,AUC_nogo,p_AUC,p_AUC_go,p_AUC_nogo,CR_ratio]=plot_AUC(act_all,CR,trial,frame);
                                        saveas(h,[savepath_figure4,figname],'png');close(h)
                                        save([savepath_figure4,figname,'.mat'],'AUC','AUC_go','AUC_nogo','p_AUC','p_AUC_go','p_AUC_nogo','CR_ratio', '-v7.3');
                                    end
                                end
                                close all;
                            end
                        end
                    else
                        warning(['Error in' savepathh]);continue;
                    end
                end
                else
                    warning(['No mat in ', savepath,'/clust_results.mat']);
                end
            end
        else
            warning(['No brain_region_related_statistic.mat in'  path_fishi])
        end
    end
end

%% 转存
clc;clear all;
seg_batchi={'Learner','Unpair','Non-learner','Faded-learner'};
load('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\Path.mat');
destination_path=checkpath(fullfile('\\10.10.42.178\DuLabSMB\13_Team Work_Update\Learning','Cluster_results_202211'));
for batchi=1:length(Path)
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi};nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
        savepath=checkpath(fullfile(path_fishi,'clust_durCond'));
        destination=checkpath(fullfile(destination_path,seg_batchi{batchi},nn));
        if exist(savepath,'dir')
            [status,msg]=copyfile(savepath,destination);
            if status==0
                disp(msg);
            end
        else
            warning([savepath ' not exist'])
        end
    end
end
%% 转存数据
clc;clear all;
seg_batchi={'Learner','Unpair','Non-learner','Faded-learner'};
mat_seg={'activities_aft_process.mat','/behavior/behav.mat','para.mat','env.mat','vol_env_spatialloc_warped_SyN_add_brainregion.csv','act_spon.mat','/behavior/behav_spon.mat'};
load('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\Path.mat');
destination_path=checkpath('J:\4.合作_20191022\数据\学习预测_ZZX');
for batchi=3:length(Path)
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi};nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
        for filei=1:length(mat_seg)
            destination=checkpath(fullfile(destination_path,seg_batchi{batchi},nn));
            savepath=fullfile(path_fishi,mat_seg{filei});
            if exist(savepath,'file')
                [status,msg]=copyfile(savepath,destination);
                if status==0
                    warning([savepath, msg]);
                end
            else
                warning([savepath ' ----------not exist'])
            end
        end
    end
end

%% 统计过程中3种类型type
%判断有无CR
%% 
outpath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
supervolxel=[];AUC_all=[];AUC_goall=[];AUC_nogoall=[];CR_all=[];fraction_all=[];AUC_all_seq=[];AUC_go_all_seq=[];AUC_nogo_all_seq=[];CR_all_seq=[];
for typei=1:length(XXX)
    path_fishi=Path{1}{double(XXX(typei,1))};nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
    Region=XXX(typei,2);
    clusteri=double(XXX(typei,3));
    
    savepath=checkpath(fullfile(path_fishi,'clust_durCond',Region))
    savepathh=fullfile(savepath,'k-means');savepath_figure=checkpath(fullfile(savepathh,'/'));
    savepath_figure4=checkpath(fullfile(savepathh,'AUC for each cluster'));
    %load cIX
    load(strcat(savepath_figure,'Fraction of clustertypes','.mat'));
    load(strcat(savepath,'/clust_results.mat'));load(strcat(savepathh,'/clust_results_s.mat'));
    cIX=cIX_s;gIX=gIX_s;
    ind=cIX(find(gIX==clusteri));
    supervolxel=[supervolxel;supervolxel_regioni(ind,:)];
    fraction_all(:,typei)=fraction(clusteri);
    %load AUC
    load(strcat(savepath_figure4,['Cluster ', num2str(clusteri)],'.mat'));
    AUC_all(:,typei)=mean(AUC,2,'omitnan');AUC=normalize(AUC,1);AUC_all_seq=[AUC_all_seq,AUC];
    AUC_goall(:,typei)=mean(AUC_go,2,'omitnan');AUC_go=normalize(AUC_go,1);AUC_go_all_seq=[AUC_go_all_seq,AUC_go];
    AUC_nogoall(:,typei)=mean(AUC_nogo,2,'omitnan');AUC_nogo=normalize(AUC_nogo,1);AUC_nogo_all_seq=[AUC_nogo_all_seq,AUC_nogo];
    CR_all(:,typei)=CR_ratio;CR_all_seq=[CR_all_seq,repmat(CR_ratio,1,size(AUC,2))];
end
AUC_all=normalize(AUC_all,1);AUC_goall=normalize(AUC_goall,1);AUC_nogoall=normalize(AUC_nogoall,1);%CR_all=normalize(CR_all,1);
AUC_all_seq_r=reshape(AUC_all_seq,1,[]);AUC_go_all_seq_r=reshape(AUC_go_all_seq,1,[]);AUC_nogo_all_seq_r=reshape(AUC_nogo_all_seq,1,[]);
CR_all_seq_r=reshape(CR_all_seq,1,[]);
clr=hsv(length(XXX));

figure,
plot(AUC_all,CR_all)

h=figure('position',[1026,735,698,243]);t=tiledlayout(h,1,3);t.TileSpacing = 'compact';t.Padding = 'compact';
for jj=1:3
    switch jj
        case 1
            AUC=AUC_all;seg='All';
        case 2
            AUC=AUC_goall;seg='Go trials';
        case 3
            AUC=AUC_nogoall;seg='No-Go trials';
    end
    ax1=nexttile(t,[1 1]);   
    for ii=1:length(XXX)
        scatter(AUC(:,ii),CR_all(:,ii),16,clr(ii,:),'filled');hold on;
    end
    [p,S]= polyfit(AUC,CR_all,1);
    [y_fit,delta] = polyval(p,AUC(:,ii),S);
    plot(AUC(:,ii),y_fit,'-','linewidth',1.5,'color','b');hold on;ylim([0 1]);xlim([-2 2])
    legend off;
    xlabel('AUC');ylabel('CR ratio');grid on;set(gca,'fontsize',14,'FontWeight','bold','linewidth',2);title(seg)
end
h=figure('position',[1026,735,698,243]);t=tiledlayout(h,1,3);t.TileSpacing = 'compact';t.Padding = 'compact';
for jj=1:3
    switch jj
        case 1
            AUC=AUC_all;seg='All';
        case 2
            AUC=AUC_goall;seg='Go trials';
        case 3
            AUC=AUC_nogoall;seg='No-Go trials';
    end
    ax1=nexttile(t,[1 1]);
    boxplot(CR_all,AUC);hold on;
end

load(fullfile('H:\3.Juvenile reference brain\registration to templete\脑区分割\','brain_region_related_statistic.mat'));
region_id=[26,29];Region=strcat(Label{region_id(1)},'-',Label{region_id(2)});
neuron_id_in_regioni=find(strcmp(brain_region_id(:,2),Label{region_id(1)}) | strcmp(brain_region_id(:,2),Label{region_id(2)}));
h=figure('position',[313,74,1099,719]);
t=tiledlayout(h,3,2);t.TileSpacing = 'compact';t.Padding = 'compact';
ax1=nexttile(t,[2 1]);
scatter3(temp_supervoxel(neuron_id_in_regioni,2),temp_supervoxel(neuron_id_in_regioni,1),temp_supervoxel(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
scatter3(supervolxel(:,2),supervolxel(:,1),supervolxel(:,3),14,'r','filled');
axis equal;view([60 40]);xlabel('X');ylabel('Y');zlabel('Z');grid on;
set(gca,'ZDir','reverse','FontWeight','bold','linewidth',2);box on;

ax1=nexttile(t,[2 1]);
scatter(temp_supervoxel(neuron_id_in_regioni,2),temp_supervoxel(neuron_id_in_regioni,1),14,[0.5 0.5 0.5],'filled');hold on;
scatter(supervolxel(:,2),supervolxel(:,1),14,'r','filled');
axis equal;;xlabel('X');ylabel('Y');grid off;
a=get(gca,'YLim');set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'ytick',[(a(2)-a(1))/2+a(1) a(2)],'yticklabel',{'Medial','Lateral'});
set(gca,'FontWeight','bold','linewidth',2);box on;view([90 90])

ax1=nexttile(t,[1 1]);
scatter(temp_supervoxel(neuron_id_in_regioni,2),temp_supervoxel(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
scatter(supervolxel(:,2),supervolxel(:,3),14,'r','filled');
axis equal;view([180 90]);xlabel('X');ylabel('Z');grid on;
legend({['n=',num2str(size(neuron_id_in_regioni,1))],['n=',num2str(length(supervolxel))]},'Location','northeast','Fontsize',10)
set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'XDir','reverse','ytick',[get(gca,'YLim')],'yticklabel',{'Dorsal','Ventral'});
set(gca,'FontWeight','bold','linewidth',2);box on;

figure,
hh=bar([1:length(XXX)]+0.5,fraction_all,'FaceColor','flat');
text([1:length(XXX)]+0.5,fraction_all,num2str(fraction_all'),'HorizontalAlignment','center','VerticalAlignment','bottom');
ylabel('Neurons');set(gca,'fontsize',14,'FontWeight','bold','linewidth',2,'xtick',1:length(XXX));box off;

im=double(temp_env.vol(:,:,1:end-1));rescalegd(im2double(temp_env.vol), [1/10000 1/10000]);
showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
for ii=1
    clr=[1 0 0];kk=1;
    zi=unique(supervolxel(:,3));
    for zz=zi'
        pti =  id(find(env.supervoxel(id,3)==zz));
        slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [supervolxel(pti, 1), supervolxel(pti, 2), floor(supervolxel(pti, 5))+3],'color',clr(kk,:));
        showspv(:,:,:,zz)=slice;
    end
    kk=kk+1;
end
seqwrite(showspv, fullfile(outpath,'emotional neurons'));
close all;