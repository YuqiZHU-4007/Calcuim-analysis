clc;clear all;
res=[0.66,0.66,10];radium=floor(30);

load('C:\run in other tap\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('C:\run in other tap\env.mat');
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
%id=find(temp_supervoxel(:,3)==(30-1)*res(3));figure,scatter(temp_supervoxel(id,1),temp_supervoxel(id,2))
mat_name=[];
Path{1}={%['H:\1.Test US\2.Tail free——Data from 117\20210402\fish1\',mat_name ],
    ['H:\1.Test US\2.Tail free——Data from 117\20210430\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free——Data from 117\20210709\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free——Data from 117\20210805\fish1\',mat_name],...
    ['E:\A_Data_lightsheet\Data_huc\20190514\fish3\',mat_name],...
    ['E:\A_Data_lightsheet\Data_huc\20190604\fish4\',mat_name],...
    ['E:\A_Data_lightsheet\Data_huc\20190609\fish3\',mat_name]
    };%learner
Path{2}={%['H:\1.Test US\2.Tail free——Data from 117\20210514\fish2\',mat_name ],...
    %['H:\1.Test US\2.Tail free——Data from 117\20210617\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free——Data from 117\20210702\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free——Data from 117\20210802\fish2\',mat_name],...
    ['E:\A_Data_lightsheet\Data_huc\20190612\fish1\',mat_name ],...
    ['E:\A_Data_lightsheet\Data_huc\20190612\fish2\',mat_name ]
    };%unpair
Path{3}={['H:\1.Test US\2.Tail free——Data from 117\20210510\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free——Data from 117\20210529\fish1\',mat_name],...
    ['G:\data_huc_non_learner\20190604\fish2\',mat_name ],...
    ['G:\data_huc_non_learner\20190607\fish3\',mat_name ]
    };%non-learner
%Path{4}={['H:\1.Test US\2.Tail free——Data from 117\20210719\fish1\',mat_name ],...
    %['H:\1.Test US\2.Tail free——Data from 117\20210717\fish1\',mat_name]};%fadd-learner
seg_batchi={'learner','unpair','non-learner','faded-learner'}; Fraction_type={};
seg1={'ON','OFF','UNRESPONSIVE'};
seg2={'Increased','Decreased','Emerged','shift to other','Increased2'};clr=jet(length(seg1)*length(seg2));
%% 按学习过程中CR聚类
addpath(genpath('C:\run in other tap\code\'));
savepathh=checkpath('C:\run in other tap\cluter_durCond\');
%单条鱼
load('C:\run in other tap\para.mat')
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
%% 
%% pre cluster 
n_pc=5;n_cluster=3;
for dr_type=2
    for clut_type=1
        savepath=checkpath(fullfile(savepathh,['n_pc ',num2str(n_pc),'--n_cluster ',num2str(n_cluster)],['dr_type ',num2str(dr_type),'--clut_type ',num2str(clut_type)]));
        seg=[['n_pc ',num2str(n_pc),'--n_cluster ',num2str(n_cluster)],['--dr_type ',num2str(dr_type),'--clut_type ',num2str(clut_type)]]
        for zz=[1 3]
            for  hh=[1 2 4]
                if exist([savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'file')
                    disp([savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'])
                    load( [savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'] );
                else
                    continue;
                end
                % figure,imagesc(M_cluster',[0,8]);colormap('hot');colorbar;title([seg '--' seg1{zz},seg2{hh}],'Interpreter','none')
                %cluster ****************************************************
                diary(fullfile(savepath,[seg1{zz},seg2{hh},'-logfile.txt']))
                [number_cluster,gIX,M_cluster_dr,Y_cluster,explained_all,n_cluster]=myclusering(M_cluster,dr_type,clut_type);
                diary off;
                save([savepath,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'M_cluster','act_all_acq_block','supervoxel_i','M_cluster_fishi',...
                   'M_cluster_dr','Y_cluster','explained_all','-v7.3');
            end
        end
    end
end
%% final cluster
Label3={'L-','C-','NL-'};
number_cluster=50;n_pc=nan;n_cluster=nan;clut_type=2;dr_type=2;
savepath=checkpath(fullfile(savepathh,['n_pc ',num2str(n_pc),'--n_cluster ',num2str(n_cluster),'--dr_type ',num2str(dr_type),'--clut_type ',num2str(clut_type)]))
seg=[['n_pc ',num2str(n_pc),'--n_cluster ',num2str(n_cluster)],['--dr_type ',num2str(dr_type),'--clut_type ',num2str(clut_type)]]
Fraction_type_durcond={};Fraction_in_region_type_durcond={}; Fraction_in_region_type_in_cluster_durcond={};
Fraction_type_durcond_fishi={}; Fraction_in_region_type_durcond_fishi={}; Fraction_in_region_type_in_cluster_durcond_fishi={};Fraction_type_durcond_fishi_typei={};Fraction_type_durcond_fishi_typei_norm={};
for zz=[1 3]
    for  hh=[1 2 4]
        if exist([savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'file')
            disp([savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'])
            load( [savepathh,'cluster_type_',seg1{zz},seg2{hh},'.mat'] );
        else
            continue;
        end
        %final cluster    
         %M_cluster_dr = tsne(M_cluster','Algorithm','exact','Distance','euclidean','Perplexity',ceil(size(M_cluster,2)/100));
        [gIX,~]=myclusering_final(M_cluster,M_cluster_dr,clut_type,number_cluster,n_pc);
        cIX=1:size(M_cluster,2);clrmap = GetColormap('hsv_new',max(gIX));
        % t-sne结果:scatter of M_cluster_dr
        figure, gscatter(M_cluster_dr(:,1),M_cluster_dr(:,2),gIX);title([seg '--' seg1{zz},seg2{hh}],'Interpreter','none');
        xlabel('Tsne-1');ylabel('Tsne-2');set(gca,'fontsize',16);legend off
        %imagesc
        figure,imagesc(M_cluster',[-2,2]);colormap('jet');colorbar;
        [~,a]=sort(gIX);b=gIX(a);
        figure,imagesc(M_cluster(:,a)',[-2,2]);colormap('jet');colorbar;title([seg '--' seg1{zz},seg2{hh}],'Interpreter','none')
        [h]=pushbutton_popupplot_Callback(M_cluster',cIX,gIX,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[savepath ,'\',seg1{zz},seg2{hh},'_popupplot'],'fig');
        %in fishi
        for ii=unique(M_cluster_fishi(2,:))
            a=find(strcmp(M_cluster_fishi(2,:),ii));
            [h]=pushbutton_popupplot_Callback(M_cluster(:,a)',1:length(a),gIX(a),clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
            text(4028,7268,ii)
            saveas(h,strcat(checkpath([savepath ,'\',seg1{zz},seg2{hh},'\','popupplot','\']),Label3{str2double(unique(M_cluster_fishi(1,a)))},ii),'png');close(h)
        end
        %plot & mapback
        for ii=unique(M_cluster_fishi(2,:))
            a=find(strcmp(M_cluster_fishi(2,:),ii));
            path=Path{str2num(unique(M_cluster_fishi(1,a)))};aa=string(path)';aa=strrep(aa,'\',''); fish_i= find(contains(aa,ii));
            if ~isempty(a);p=path{fish_i};
                env=load(fullfile(p,'env.mat'));
            else
                env.env.supervoxel=[];
            end
            supervoxel=[env.env.supervoxel(:,2) env.env.supervoxel(:,1) env.env.supervoxel(:,3)];
            h=DrawTiledPics_zyq_20190530(str2double(M_cluster_fishi(3,a)),gIX(a),[1:size(supervoxel,1)],supervoxel,env.env.vol,clrmap);
             saveas(h,strcat(checkpath([savepath ,'\',seg1{zz},seg2{hh},'\','mapback_fishi','\']),Label3{str2double(unique(M_cluster_fishi(1,a)))},ii),'png');close(h)
            text(10,10,ii,'color','w')
        end
        for ii=unique(M_cluster_fishi(2,:))
            a=find(strcmp(M_cluster_fishi(2,:),ii));
            path=Path{str2num(unique(M_cluster_fishi(1,a)))};aa=string(path)';aa=strrep(aa,'\',''); fish_i= find(contains(aa,ii));
            if ~isempty(a);p=path{fish_i};
                env=load(fullfile(p,'env.mat'));
            else
                env.env.supervoxel=[];
            end
            supervoxel=[env.env.supervoxel(:,2) env.env.supervoxel(:,1) env.env.supervoxel(:,3)];
            h=figure('position',[7,154,1910,824]);
            t=tiledlayout(ceil(length(unique(gIX(a)))/10),10);t.TileSpacing = 'none';t.Padding = 'none';
            for jj=unique(gIX(a))'
                aa=a(find(gIX(a)==jj));
                nexttile,scatter3(supervoxel(:,1)*res(1),supervoxel(:,2)*res(2),(supervoxel(:,3)-1)*res(3),4,[0.5 0.5 0.5]);hold on
                aaa=str2double(M_cluster_fishi(3,aa));
                scatter3(supervoxel(aaa,1)*res(1),supervoxel(aaa,2)*res(2),(supervoxel(aaa,3)-1)*res(3),12,clrmap(gIX(aa),:),'filled');hold on;
                axis equal;set(gca,'xtick',[],'ytick',[]);view([-90 90]);
                text(gca,400,1200,num2str(length(aa)))
                title(num2str(jj));
                if jj==min(unique(gIX(a))')
                    text(1100,1250,ii);
                end
            end 
             saveas(h,strcat(checkpath([savepath ,'\',seg1{zz},seg2{hh},'\','mapback_scatter_fishi','\']),Label3{str2double(unique(M_cluster_fishi(1,a)))},ii),'png');close(h)
        end
        for ii=unique(M_cluster_fishi(2,:))
            a=find(strcmp(M_cluster_fishi(2,:),ii));
            path=Path{str2num(unique(M_cluster_fishi(1,a)))};aa=string(path)';aa=strrep(aa,'\',''); fish_i= find(contains(aa,ii));
            if ~isempty(a);p=path{fish_i};
                env=load(fullfile(p,'env.mat'));
            else
                env.env.supervoxel=[];
            end
            supervoxel=[env.env.supervoxel(:,2) env.env.supervoxel(:,1) env.env.supervoxel(:,3)];
            for jj=unique(gIX(a))'
                aa=a(find(gIX(a)==jj));
                aaa=str2double(M_cluster_fishi(3,aa));
                h=DrawTiledPics_zyq_20190530(aaa,gIX(aa),[1:size(supervoxel,1)],supervoxel,env.env.vol,clrmap);
                text(10,10,strcat(ii,'-',num2str(jj)),'color','w');
                saveas(h,strcat(checkpath(strcat(savepath ,'\',seg1{zz},seg2{hh},'\','mapback_scatter_fishi_typei','\',Label3{str2double(unique(M_cluster_fishi(1,a)))},ii,'\')),num2str(jj)),'png');close(h)
            end
        end
        %1
        [h,ratio]=plot_test_2(act_all_acq_block,cIX,gIX,frame,clrmap,4,[-2 2],true,true);saveas(h,[savepath ,'\',seg1{zz},seg2{hh},'_trace'],'fig');
        %2
        figure('position',[7,154,1910,824]),t=tiledlayout(1,3);t.TileSpacing = 'compact';t.Padding = 'compact';
        for batchi=1:3
            ind_i=find(str2double(M_cluster_fishi(1,:))==batchi);
            supervoxel=supervoxel_i(ind_i,:);gIX_i=gIX(ind_i,:);disp([seg_batchi{batchi},'-supervoxel number:',num2str(size(supervoxel,1))])
            nexttile,scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),10,clrmap(gIX_i,1:3),'filled');
            axis equal;ylim([1 temp_env.env.height*res(1)]);xlim([1 temp_env.env.width*res(1)]);zlim([0 43]*res(3));view(180,90);grid on;%view([138 58]);
        end
        saveas(h,[savepath,'\' ,seg1{zz},seg2{hh},'_scatter mapback'],'fig');
        %2
        figure('position',[7,154,1910,824]),t=tiledlayout(1,3);t.TileSpacing = 'compact';t.Padding = 'compact';
        for batchi=1:3
            ind_i=find(str2double(M_cluster_fishi(1,:))==batchi);supervoxel=supervoxel_i(ind_i,:);gIX_i=gIX(ind_i,:);
            [~,ind_cut,ind_save,~]=cut_small_cluster(gIX_i,20);disp([seg_batchi{batchi},'-supervoxel number cut:',num2str(length(ind_save))])
            nexttile, scatter3(supervoxel(ind_save,1),supervoxel(ind_save,2),supervoxel(ind_save,3),10,clrmap(gIX_i(ind_save),1:3),'filled');
            hold on;scatter3(supervoxel(ind_cut,1),supervoxel(ind_cut,2),supervoxel(ind_cut,3),10,[0.5 0.5 0.5],'filled');
            axis equal;ylim([1 temp_env.env.height*res(1)]);xlim([1 temp_env.env.width*res(1)]);zlim([0 43]*res(3));view([138 58]);grid on;
        end
        saveas(h,[savepath,'\' ,seg1{zz},seg2{hh},'_scatter mapback_cut'],'fig');
        %脑区占比
        for batchi=1:3
            path=Path{batchi};aa=string(path)';aa=strrep(aa,'\','');
            ind_i=find(str2double(M_cluster_fishi(1,:))==batchi );supervoxel=supervoxel_i(ind_i,:);gIX_i=gIX(ind_i,:);
            Fraction_type_durcond{batchi,zz}(hh)=length(ind_i)./size(temp_supervoxel,1);
            [fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust]=get_region_fraction_temp_preprocess(gIX_i,supervoxel,reg_mask,reg_name,reg_loc,temp_supervoxel,temp_env,clrmap);
            Fraction_in_region_type_durcond{batchi,zz}(:,hh,1:max(gIX))=fraction_in_region_in_clust;
            Fraction_in_region_type_in_cluster_durcond{1}{batchi,zz}(:,hh,1:max(gIX))= num_in_region_in_clust{1};
            Fraction_in_region_type_in_cluster_durcond{2}{batchi,zz}(:,hh,1:max(gIX))= num_in_region_in_clust{2};
            Fraction_in_region_type_in_cluster_durcond{3}{batchi,zz}(:,hh,1:max(gIX))= num_in_region_in_clust{3};
            kkkk=1;
            for ii=unique(M_cluster_fishi(2,:))
                if ~isempty(strfind(ii,'2019'))
                    res2=[0.66,0.66,8];
                elseif ~isempty(strfind(ii,'2021'))
                    res2=[0.66,0.66,10];
                end
                fish_i= find(contains(aa,ii));
                ind_i= find(strcmp(ii,M_cluster_fishi(2,:)) & str2double(M_cluster_fishi(1,:))==batchi);
                supervoxel=supervoxel_i(ind_i,:);gIX_i=gIX(ind_i,:);
                if ~isempty(fish_i);p=path{fish_i};
                    env=load(fullfile(p,'env.mat')); 
                else
                    env.env.supervoxel=[];
                end
                Fraction_type_durcond_fishi{batchi,zz}(hh,kkkk)=length(ind_i)./size(env.env.supervoxel,1);
                for jj=unique(gIX)'
                    ind_i= find(strcmp(ii,M_cluster_fishi(2,:)) & str2double(M_cluster_fishi(1,:))==batchi & gIX'==jj);
                    ind_fishi= find(strcmp(ii,M_cluster_fishi(2,:)) & str2double(M_cluster_fishi(1,:))==batchi);
                    Fraction_type_durcond_fishi_typei{batchi,zz}(hh,kkkk,jj)=length(ind_i)./size(env.env.supervoxel,1);
                    Fraction_type_durcond_fishi_typei_norm{batchi,zz}(hh,kkkk,jj)=length(ind_i)./length(ind_fishi);
                end
 
                [fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust]=get_region_fraction_temp_preprocess(gIX_i,supervoxel,reg_mask,reg_name,reg_loc,temp_supervoxel,temp_env,clrmap);
                Fraction_in_region_type_durcond_fishi{batchi,zz}(:,hh,kkkk,1:size(fraction_in_region_in_clust,2))=fraction_in_region_in_clust;
                Fraction_in_region_type_in_cluster_durcond_fishi{1}{batchi,zz}(:,hh,kkkk,1:size(fraction_in_region_in_clust,2))= num_in_region_in_clust{1};
                Fraction_in_region_type_in_cluster_durcond_fishi{2}{batchi,zz}(:,hh,kkkk,1:size(fraction_in_region_in_clust,2))= num_in_region_in_clust{2};
                Fraction_in_region_type_in_cluster_durcond_fishi{3}{batchi,zz}(:,hh,kkkk,1:size(fraction_in_region_in_clust,2))= num_in_region_in_clust{3};
                kkkk=kkkk+1;
            end
        end  
        %save
        save([savepath,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'M_cluster','act_all_acq_block','supervoxel_i','M_cluster_fishi',...
            'number_cluster','cIX' ,'gIX','M_cluster_dr','Y_cluster','explained_all','clrmap','n_pc',...
            'Fraction_type_durcond','Fraction_in_region_type_durcond', 'Fraction_in_region_type_in_cluster_durcond',...
            'Fraction_type_durcond_fishi', 'Fraction_in_region_type_durcond_fishi', 'Fraction_in_region_type_in_cluster_durcond_fishi','Fraction_type_durcond_fishi_typei','Fraction_type_durcond_fishi_typei_norm','-v7.3');
        
    end
end
%% 手动统计
%脑区各类占比 /脑区内总neuron
Lable={'L OB','R OB','L P','R P',...
    'L Hb','R Hb','L TeO','R TeO','L Np','R Np','L Otr','R Otr','TL','PO',...
    'Pr','PT','Th','rH','iH','cH','mPT','T',...
    'L TS','R TS','Va','CCe','L gT','R gT','R1','R2'};
%Label2={'Activation-Increase','Activation-Decrease','Activation-Shift','Inhibition-Increase','Inhibition-Decrease','Inhibition-Shift','Activation-Emerged','Inhibition-Emerged'};
Label3={'Learner','Control','Non-learner'};
%zz=1;hh=1;
k=2;
for plot_type=2
    for type_i=unique(gIX)'
        aa=[];
        for batchi=1:3
            a=(squeeze( Fraction_in_region_type_durcond_fishi{batchi,zz}(:,hh,:,type_i))*100);
            aa=[aa,mean(a,2,'omitnan')];
        end
        if mod(type_i,ceil(max(gIX)/k))==1
            figure('position',[0,0,758,1000]),t=tiledlayout(ceil(max(gIX)/k),1);t.TileSpacing = 'compact';t.Padding = 'compact';           
        end
        nexttile,b=bar(aa,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);title(num2str(type_i));
        if mod(type_i,ceil(max(gIX)/k))==0 %type_i>=max(gIX)-k
            %set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',15);          
            ylabel('Fraction');legend({'L','C','NL'});set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',10);
            legend({'Learner','Control','Non-learner'});
        elseif type_i==1
            %ylabel('Fraction'); set(gca,'xtick',[1:length(Lable)],'xticklabels',{},'fontsize',10);
        else
            set(gca,'xtick',[1,length(Lable)],'xticklabels',{},'fontsize',10);
            %yaxis_cut(aa)
        end
    end
end
Label2=mat2cell([1:max(gIX)],1,ones(max(gIX),1));
for plot_type=2
    aa=[];bb=[];
    for type_i=unique(gIX)'
        for batchi=1:3
            a=squeeze(Fraction_type_durcond_fishi_typei{batchi,zz}(hh,:,type_i)*100);
            b=squeeze(Fraction_type_durcond_fishi_typei_norm{batchi,zz}(hh,:,type_i));
            aa(batchi,type_i,:)=a;
            bb(batchi,type_i,:)=b;
        end
    end
    figure('position',[46,214,1400,479]),
    b=bar(mean(aa,3,'omitnan')','stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label2)]);legend({'Learner','Control','Non-learner'});
    set(gca,'xtick',[1:length(Label2)],'xticklabels',Label2,'XTickLabelRotation',45,'fontsize',15);
    
    figure('position',[46,43,337,933]),
    b=bar(mean(aa,3,'omitnan'),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label3)]);kk=1;
    for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
        b(k).CData = clrmap(kk,:);
        kk=kk+1;
    end
    set(gca,'xtick',[1:length(Label3)],'xticklabels',Label3,'XTickLabelRotation',45,'fontsize',15);
    
    figure('position',[46,43,337,933]),
    b=bar(mean(aa,3,'omitnan'),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label3)]);kk=1;
    for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
        b(k).CData = clrmap(kk,:);
        kk=kk+1;
    end
    set(gca,'xtick',[1:length(Label3)],'xticklabels',Label3,'XTickLabelRotation',45,'fontsize',15);

    figure('position',[46,43,337,933]),
    b=bar(mean(bb,3,'omitnan'),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label3)]);kk=1;
    for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
        b(k).CData = clrmap(kk,:);
        kk=kk+1;
    end
    ylim([0 1]);set(gca,'xtick',[1:length(Label3)],'xticklabels',Label3,'XTickLabelRotation',45,'fontsize',15);
end
for plot_type=2
    %t=tiledlayout(3,1);t.TileSpacing = 'compact';t.Padding = 'compact';nexttile,
    for batchi=1:3
        a=squeeze(Fraction_in_region_type_durcond_fishi{batchi,zz}(:,hh,:,:)*100);
        aa=squeeze(mean(a,2,'omitnan'));
        figure('position',[0,0,900,600]),
        b=bar(aa,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);kk=1;
        for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
            b(k).CData = clrmap(kk,:);
            kk=kk+1;
        end
        title(Label3{batchi});
        set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',15);
    end
end
for plot_type=2
    for batchi=1:3
        a=squeeze(Fraction_in_region_type_in_cluster_durcond_fishi{2}{batchi,zz}(:,hh,:,:));
        b=squeeze(sum(a,3,'omitnan'));
        c=a./repmat(b,1,1,size(a,3));
         figure('position',[0,0,900,600]),
         b=bar(squeeze(mean(c,2,'omitnan')),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);kk=1;
        for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
            b(k).CData = clrmap(kk,:);
            kk=kk+1;
        end
        title(Label3{batchi});
        ylim([0 1]);set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',15);
    end
end
%% 取ind
ind=[25,48,23,19, 6, 46, 13, 22,9, 49,20, 31, 47,  16,18,15,17, 26,32, 37, 38,39,45, 44, 12,40];clrmap_ind=clrmap(ind,:);
kk=1;gIX_ind=[];gIX_new=[];
for ii=1:length(ind)
    gIX_indi=find(gIX==ind(ii));gIX_ind=[gIX_ind;gIX_indi];
    gIX_new=[gIX_new,kk*ones(1,length(gIX_indi))];kk=kk+1;
end
gIX(gIX_ind)
figure,imagesc(M_cluster(:,gIX_ind)',[0,6]);colormap('hot');colorbar;
[h]=pushbutton_popupplot_Callback(M_cluster(:,gIX_ind)',1:length(gIX_ind),gIX_new,clrmap_ind,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
for plot_type=2
    figure('position',[0,0,758,1000]),t=tiledlayout(ceil(length(ind)/2),2);t.TileSpacing = 'compact';t.Padding = 'compact';   
    for ii=1:length(ind)
        aa=[];
        for batchi=1:3
            a=(squeeze( Fraction_in_region_type_durcond_fishi{batchi,zz}(:,hh,:,ind(ii)))*100);
            aa(:,batchi)=[mean(a,2,'omitnan')];
        end
        nexttile,b=bar(aa,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);title(num2str(ind(ii)));
        if ii==length(ind)     
            ylabel('Fraction');legend({'L','C','NL'});set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',10);
            legend({'Learner','Control','Non-learner'});
        else
            set(gca,'xtick',[1,length(Lable)],'xticklabels',{},'fontsize',10);
        end
    end
end
Label2=mat2cell(ind,1,ones(length(ind),1));
for plot_type=2
    aa=[];bb=[];
    for ii=1:length(ind)
        type_i=ind(ii);
        for batchi=1:3
            a=squeeze(Fraction_type_durcond_fishi_typei{batchi,zz}(hh,:,type_i)*100);
            b=squeeze(Fraction_type_durcond_fishi_typei_norm{batchi,zz}(hh,:,type_i));
            aa(batchi,ii,1:length(a))=a;
            bb(batchi,ii,1:length(a))=b;
        end
    end
    figure('position',[46,214,1400,479]),
    b=bar(mean(aa,3,'omitnan')','stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label2)]);legend({'Learner','Control','Non-learner'});
    set(gca,'xtick',[1:length(Label2)],'xticklabels',Label2,'XTickLabelRotation',45,'fontsize',15);
    
    figure('position',[46,43,337,933]),
    b=bar(mean(aa,3,'omitnan'),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label3)]);
    for k =1:length(ind)
        b(k).CData = clrmap(ind(k),:);
    end
    set(gca,'xtick',[1:length(Label3)],'xticklabels',Label3,'XTickLabelRotation',45,'fontsize',15); 
    
    figure('position',[46,43,337,933]),
    b=bar(mean(bb,3,'omitnan'),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Label3)]);kk=1;
    for k =1:length(ind)
        b(k).CData = clrmap(ind(k),:);
    end
    ylim([0 1]);set(gca,'xtick',[1:length(Label3)],'xticklabels',Label3,'XTickLabelRotation',45,'fontsize',15);
end
for plot_type=2
    %t=tiledlayout(3,1);t.TileSpacing = 'compact';t.Padding = 'compact';nexttile,
    for batchi=1:3
        a=squeeze(Fraction_in_region_type_durcond_fishi{batchi,zz}(:,hh,:,ind)*100);
        aa=squeeze(mean(a,2,'omitnan'));
        figure('position',[0,0,900,600]),
        b=bar(aa,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);
        for k =1:length(ind)
            b(k).CData = clrmap(ind(k),:);
        end
        title(Label3{batchi});
        set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',15);
    end
end
for plot_type=2
    for batchi=1:3
        a=squeeze(Fraction_in_region_type_in_cluster_durcond_fishi{2}{batchi,zz}(:,hh,:,ind));
        b=squeeze(sum(a,3,'omitnan'));
        c=a./repmat(b,1,1,size(a,3));
         figure('position',[0,0,900,600]),
         b=bar(squeeze(mean(c,2,'omitnan')),'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);
        for k =1:length(ind)
            b(k).CData = clrmap(ind(k),:);
        end
        title(Label3{batchi});
        ylim([0 1]);set(gca,'xtick',[1:length(Lable)],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',15);
    end
end
