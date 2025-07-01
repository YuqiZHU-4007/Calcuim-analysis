function plot_area_from_excel
%写每列数据到excel方便graphpad画图；画每个脑区不同鱼不同列数据的图
clc;clear all;close all;
set(0,'defaultfigurecolor','w');
%可修改parameters
global XLrange %area_mean 写入excel第几列，新家数据一定要改，不然会覆盖之前一条鱼的数据！！！
global n %matric：提取每条鱼表中第几列area_mean
global paranum
global rownum_eachfish
global isflash
global norm %logic：是否normalize
global normmethod %char：normalize 方法
global figure_able %logic：是否画图
global figure_able_mean
global figure_able_contolCS %logic ：是否画learner部分
global figure_able_contolUS %logic ：是否画control CS部分
global figure_able_learner %logic ：是否画control US部分
global lab %logic：是否需要legend

global linewidth
global marksize
global color_line
global color_mark
global color_line_controlCS
global color_mark_controlCS
global color_line_controlUS
global color_mark_controlUS

%%%control和 learner的sheet name要相同
filepath='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328\massive\area_event';
controlpath_CS='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328\control-CS\area_event';
controlpath_US='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328\control-US\area_event';

%para.
isflash=true;
fishnum=4:5;
paranum=9;
rownum_eachfish=paranum+1;
n=4;%读取第几列参数
norm=false;
normmethod='zscore';%'range','zscore','mean'(减去均值),'first'（减去第一个block数值）
lab=true;
figure_able=true;
figure_able_mean=true;
figure_able_learner=true;
figure_able_contolCS=false;
figure_able_contolUS=false;


linewidth=2.5;
marksize=5;
color_line=[];color_mark='r';
color_line_controlCS=[0.5 0.5 0.5];color_mark_controlCS='k';
color_line_controlUS=[0.5 0 0];color_mark_controlUS='b';


[~,sheet,~]=xlsfinfo(filepath);
name_fishnumber={}; name_fishnumber_CS={};name_fishnumber_US={};
% C={};%C=[];
% C1=[];C2=[];
% for ii=1:max(size(sheet))
%     sheet{ii}
%     [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
%     [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
%     %画learner 和 写入
%     [data,~,raw_data]=xlsread(filepath,sheet{ii});
%     %normalize
%     raw_data=cell2mat(raw_data(4:end,:));
%     raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
%     raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
%     raw_data=cur_nan_col(raw_data);
%     raw_data_controlCS=cur_nan_col(raw_data_controlCS);
%     raw_data_controlUS=cur_nan_col(raw_data_controlUS);
%     %extra data
%     A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
%     [A,B]=get_X_from_raw(n,raw_data,rownum_eachfish);
%     [A_control_CS,B_control_CS]=get_X_from_raw(n,raw_data_controlCS,rownum_eachfish);
%     [A_control_US,B_control_US]=get_X_from_raw(n,raw_data_controlUS,rownum_eachfish);
%     if norm
%         B=normlize_zyq_20190320(n,B,normmethod);
%         B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
%         B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
%     end
%     %C(ii,:)=[B(1,:),B(2,:)];
%     C(ii,:)=[B(1,:);B(end,:)]'
%     C1(ii,1:size(B,2))=B(1,:);
%     C2(ii,1:size(B,2))=B(end,:);
% end


a=-0.05;b=0.2;ax=0.5;bx=20.5;C={};
linewidth=2.5;
marksize=5;
color_line=[0.31 0.31 0.31];color_mark=[0.31 0.31 0.31];
color_line_controlCS=[0.5 0.5 0.5];color_mark_controlCS='k';
color_line_controlUS=[0.5 0 0];color_mark_controlUS='b';
h1=figure;%ax1=axes(h1);
for ii=1:max(size(sheet))
    [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
    [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
    %画learner 和 写入
    [data,~,raw_data]=xlsread(filepath,sheet{ii});
    name_fishnumber=get_fishnumber(raw_data,rownum_eachfish);
    name_fishnumber_CS=get_fishnumber(raw_data_controlCS,rownum_eachfish);
    name_fishnumber_US=get_fishnumber(raw_data_controlUS,rownum_eachfish);
    %normalize
    raw_data=cell2mat(raw_data(4:end,:)); raw_dataa=[];
     for kkk=1:length(fishnum) raw_dataa=[raw_dataa, raw_data(:,(fishnum(kkk)-1)*rownum_eachfish+1:min(fishnum(kkk)*rownum_eachfish,size(raw_data,2)))]; end
     raw_data=raw_dataa;
    raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
    raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
    raw_data=cur_nan_col(raw_data);
    raw_data_controlCS=cur_nan_col(raw_data_controlCS);
    raw_data_controlUS=cur_nan_col(raw_data_controlUS);
    %extra data
    A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
    [A,B]=get_X_from_raw(n,raw_data,rownum_eachfish);
    [A_control_CS,B_control_CS]=get_X_from_raw(n,raw_data_controlCS,rownum_eachfish);
    [A_control_US,B_control_US]=get_X_from_raw(n,raw_data_controlUS,rownum_eachfish);
    if norm
        B=normlize_zyq_20190320(n,B,normmethod);
        B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
        B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
    end
    %C=B(:,1:2:size(B,2))';
    if length(n)==1
        %画图
        if figure_able
            subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii),
            plot_mean_area(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold off;
            set(gca,'linewidth',2,'fontsize',20);
        end
    elseif length(n)==2 %n为 m and sd
        %画图
        if figure_able
            subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii),
            plot_mean_area_errorbar(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold off;
            set(gca,'linewidth',2,'fontsize',20);
        end
    end%n为 m and sd
    xlim([ax bx]);
    ylim([a b]);
    %     %写入data
    %     [status,message]=xlswrite(outfilepath,A,sheet{ii},'B1');
    %     %判断写入
    %     if status==0
    %         error(message);
    %     elseif status==1
    %         disp([ sheet{ii} ' down']);
    %     end
    B=cur_nan_col(B);
    C{ii}=B';
end

h2=figure;%ax2=axes(h2);
for ii=1:max(size(sheet))
    [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
    [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
    %画learner 和 写入
    [data,~,raw_data]=xlsread(filepath,sheet{ii});
    name_fishnumber=get_fishnumber(raw_data,rownum_eachfish);
    name_fishnumber_CS=get_fishnumber(raw_data_controlCS,rownum_eachfish);
    name_fishnumber_US=get_fishnumber(raw_data_controlUS,rownum_eachfish);
    %normalize
    raw_data=cell2mat(raw_data(4:end,:));
    raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
    raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
    raw_data=cur_nan_col(raw_data);
    raw_data_controlCS=cur_nan_col(raw_data_controlCS);
    raw_data_controlUS=cur_nan_col(raw_data_controlUS);
    %extra data
    A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
    [A,B]=get_X_from_raw(n,raw_data,rownum_eachfish);
    if isflash
        B(:,1:length(n))=[];
    end
    [A_control_CS,B_control_CS]=get_X_from_raw(n,raw_data_controlCS,rownum_eachfish);
    [A_control_US,B_control_US]=get_X_from_raw(n,raw_data_controlUS,rownum_eachfish);
    if norm
        B=normlize_zyq_20190320(n,B,normmethod);
        B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
        B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
    end
    if figure_able_mean
        if isflash && ~isempty(intersect(n,[1,2,4,5,6,7])) && (size(B,2)==size(name_fishnumber,2) || size(B,2)==size(name_fishnumber,2)*2)
            B(:,1:length(n))=[];
        end
        meanB=[];meanB_CS=[];meanB_US=[];
        if length(n)==1
            meanB(:,1)=mean(B,2,'omitnan');meanB(:,2)=std(B,[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS,2,'omitnan');meanB_CS(:,2)=std(B_control_CS,[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US,2,'omitnan');meanB_US(:,2)=std(B_control_US,[],2,'omitnan');
        elseif length(n)==2
            meanB(:,1)=mean(B(:,1:2:size(B,2)),2,'omitnan');meanB(:,2)=std(B(:,1:2:size(B,2)),[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS(:,1:2:size(B_control_CS,2)),2,'omitnan');meanB_CS(:,2)=std(B_control_CS(:,1:2:size(B_control_CS,2)),[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US(:,1:2:size(B_control_US,2)),2,'omitnan');meanB_US(:,2)=std(B_control_US(:,1:2:size(B_control_US,2)),[],2,'omitnan');
        end
        subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii);
        plot_mean_area_errorbar(meanB,meanB_CS,meanB_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);
        set(gca,'linewidth',2,'fontsize',20);
        text_n(n,B,B_control_CS,B_control_US,figure_able,figure_able_contolCS,figure_able_contolUS);
         xlim([ax bx]); ylim([a b]);
    end
end

h3=figure;%hab，acq，test连在一起
n1=8:9;%hab_tst
n2=4:5;%acq
% if  ~isempty(intersect(n1,[1,2,4,5,6,7]))
%     X=[ 'Acq.-Bloc.1'; 'Acq.-Bloc.2'; 'Acq.-Bloc.3'; 'Acq.-Bloc.4' ;'Acq.-Bloc.5'; ...
%         'Acq.-Bloc.6'; 'Acq.-Bloc.7'; 'Acq.-Bloc.8'; 'Acq.-Bloc.9'; 'Acq.-Bloc.10']';
% else
%     X=[['Before'+newline+'learning'] ['After'+newline+ 'learning']]';
% end
for ii=1:max(size(sheet))
    [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
    [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
    %画learner 和 写入
    [data,~,raw_data]=xlsread(filepath,sheet{ii});
    name_fishnumber=get_fishnumber(raw_data,rownum_eachfish);
    name_fishnumber_CS=get_fishnumber(raw_data_controlCS,rownum_eachfish);
    name_fishnumber_US=get_fishnumber(raw_data_controlUS,rownum_eachfish);
    %normalize
    raw_data=cell2mat(raw_data(4:end,:));
    raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
    raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
    raw_data=cur_nan_col(raw_data);
    raw_data_controlCS=cur_nan_col(raw_data_controlCS);
    raw_data_controlUS=cur_nan_col(raw_data_controlUS);
    %extra data
    A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
    [A1,B1]=get_X_from_raw(n1,raw_data,rownum_eachfish);
    [A_control_CS1,B_control_CS1]=get_X_from_raw(n1,raw_data_controlCS,rownum_eachfish);
    [A_control_US1,B_control_US1]=get_X_from_raw(n1,raw_data_controlUS,rownum_eachfish);
    [A2,B2]=get_X_from_raw(n2,raw_data,rownum_eachfish);
    [A_control_CS2,B_control_CS2]=get_X_from_raw(n2,raw_data_controlCS,rownum_eachfish);
    [A_control_US2,B_control_US2]=get_X_from_raw(n2,raw_data_controlUS,rownum_eachfish);
    A=get_integrate_X1_X2(A1,A2);B=get_integrate_X1_X2(B1,B2);
    A_control_CS=get_integrate_X1_X2(A_control_CS1,A_control_CS2);B_control_CS=get_integrate_X1_X2(B_control_CS1,B_control_CS2);
    A_control_US=get_integrate_X1_X2(A_control_US1,A_control_US2);B_control_US=get_integrate_X1_X2(B_control_US1,B_control_US2);
    if norm
        B=normlize_zyq_20190320(n,B,normmethod);
        B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
        B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
    end
    if figure_able %画图
        AX=subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii);
        Y=[B(:);B_control_CS(:);B_control_US(:)];
        line([2 2],[a b],'linewidth',2,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        line([size(B2,1)+1 size(B2,1)+1],[a b],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        %         patch([2 size(B2,1)+1 size(B2,1)+1 2  ],[AX.YLim(1)   AX.YLim(1)  AX.YLim(2) AX.YLim(2)],[0.5 0.5 0.5],...
        %             'edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.3);hold on;
        if length(n1)==1
            plot_mean_area(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);
            hold off;
        elseif length(n1)==2 %n为 m and sd
            plot_mean_area_errorbar(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold off;
        end
        %         line([2 2],[AX.YLim(1) AX.YLim(2)],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        %         line([size(B2,1)+1 size(B2,1)+1],[AX.YLim(1) AX.YLim(2)],'linewidth',1.5,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        AX.XTickLabel =[[]; 'Hab.';string([1:size(B2,1)]');'Tst.'];
        set(gca,'linewidth',2,'fontsize',20);
        ylim([a b]);
        xlim([ax bx]);
    end
end %n为 m and sd


%add mean
linewidth=2.5;
marksize=4;
color_line=[0.5 0.5 0.5];color_mark='k';
color_line_controlCS=[0.5 0.5 0.5];color_mark_controlCS='k';
color_line_controlUS=[0.5 0 0];color_mark_controlUS='b';
a=-0.05;b=0.405;
h4=figure;%ax1=axes(h1);
for ii=1:max(size(sheet))
    [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
    [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
    %画learner 和 写入
    [data,~,raw_data]=xlsread(filepath,sheet{ii});
    name_fishnumber=get_fishnumber(raw_data,rownum_eachfish);
    name_fishnumber_CS=get_fishnumber(raw_data_controlCS,rownum_eachfish);
    name_fishnumber_US=get_fishnumber(raw_data_controlUS,rownum_eachfish);
    %normalize
    raw_data=cell2mat(raw_data(4:end,:));
    raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
    raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
    raw_data=cur_nan_col(raw_data);
    raw_data_controlCS=cur_nan_col(raw_data_controlCS);
    raw_data_controlUS=cur_nan_col(raw_data_controlUS);
    %extra data
    A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
    [A,B]=get_X_from_raw(n,raw_data,rownum_eachfish);
    [A_control_CS,B_control_CS]=get_X_from_raw(n,raw_data_controlCS,rownum_eachfish);
    [A_control_US,B_control_US]=get_X_from_raw(n,raw_data_controlUS,rownum_eachfish);
    if norm
        B=normlize_zyq_20190320(n,B,normmethod);
        B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
        B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
    end
    if length(n)==1
        %画图
        if figure_able
            subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii),
            plot_mean_area(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold on;
            if isflash && ~isempty(intersect(n,[1,2,4,5,6,7])) && (size(B,2)==size(name_fishnumber,2) || size(B,2)==size(name_fishnumber,2)*2)
                B(:,1:length(n))=[];
            end
            meanB=[];meanB_CS=[];meanB_US=[];
            meanB(:,1)=mean(B(:,1:2:size(B,2)),2,'omitnan');meanB(:,2)=std(B(:,1:2:size(B,2)),[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS(:,1:2:size(B_control_CS,2)),2,'omitnan');meanB_CS(:,2)=std(B_control_CS(:,1:2:size(B_control_CS,2)),[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US(:,1:2:size(B_control_US,2)),2,'omitnan');meanB_US(:,2)=std(B_control_US(:,1:2:size(B_control_US,2)),[],2,'omitnan');
            plot_mean_area_errorbar(meanB,meanB_CS,meanB_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);
            text_n(n,B,B_control_CS,B_control_US,figure_able,figure_able_contolCS,figure_able_contolUS);
            set(gca,'linewidth',2,'fontsize',20);
            ylim([-0.05 0.2]);
        end
    elseif length(n)==2 %n为 m and sd
        %画图
        if figure_able
            subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii),
            [p1,p2,p3]=plot_mean_area(B(:,1:2:size(B,2)),B_control_CS(:,1:2:size(B_control_CS,2)),B_control_US(:,1:2:size(B_control_US,2)),...
                sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold on;
            for zz=1:length(p1) p1(zz).Color=color_line; end
            if isflash && ~isempty(intersect(n,[1,2,4,5,6,7])) && (size(B,2)==size(name_fishnumber,2) || size(B,2)==size(name_fishnumber,2)*2)
                B(:,1:length(n))=[];
            end
            meanB=[];meanB_CS=[];meanB_US=[];
            meanB(:,1)=mean(B(:,1:2:size(B,2)),2,'omitnan');meanB(:,2)=std(B(:,1:2:size(B,2)),[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS(:,1:2:size(B_control_CS,2)),2,'omitnan');meanB_CS(:,2)=std(B_control_CS(:,1:2:size(B_control_CS,2)),[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US(:,1:2:size(B_control_US,2)),2,'omitnan');meanB_US(:,2)=std(B_control_US(:,1:2:size(B_control_US,2)),[],2,'omitnan');
            [p1,p2,p3]=plot_mean_area_errorbar(meanB,meanB_CS,meanB_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold on;
            for zz=1:length(p1) p1(zz).Color='r';p1(zz).LineWidth=3;p1(zz).MarkerEdgeColor='r';p1(zz).MarkerFaceColor='r'; end
            text_n(n,B,B_control_CS,B_control_US,figure_able,figure_able_contolCS,figure_able_contolUS);
            set(gca,'linewidth',2,'fontsize',20);
            ylim([a b]);
        end
    end %n为 m and sd
    %     %写入data
    %     [status,message]=xlswrite(outfilepath,A,sheet{ii},'B1');
    %     %判断写入
    %     if status==0
    %         error(message);
    %     elseif status==1
    %         disp([ sheet{ii} ' down']);
    %     end
end

h5=figure;%hab，acq，test连在一起
n1=8:9;%hab_tst
n2=4:5;%acq
% if  ~isempty(intersect(n1,[1,2,4,5,6,7]))
%     X=[ 'Acq.-Bloc.1'; 'Acq.-Bloc.2'; 'Acq.-Bloc.3'; 'Acq.-Bloc.4' ;'Acq.-Bloc.5'; ...
%         'Acq.-Bloc.6'; 'Acq.-Bloc.7'; 'Acq.-Bloc.8'; 'Acq.-Bloc.9'; 'Acq.-Bloc.10']';
% else
%     X=[['Before'+newline+'learning'] ['After'+newline+ 'learning']]';
% end
for ii=1:max(size(sheet))
    [~,~,raw_data_controlCS]=xlsread(controlpath_CS,sheet{ii});
    [~,~,raw_data_controlUS]=xlsread(controlpath_US,sheet{ii});
    %画learner 和 写入
    [data,~,raw_data]=xlsread(filepath,sheet{ii});
    name_fishnumber=get_fishnumber(raw_data,rownum_eachfish);
    name_fishnumber_CS=get_fishnumber(raw_data_controlCS,rownum_eachfish);
    name_fishnumber_US=get_fishnumber(raw_data_controlUS,rownum_eachfish);
    %normalize
    raw_data=cell2mat(raw_data(4:end,:));
    raw_data_controlCS=cell2mat(raw_data_controlCS(4:end,:));
    raw_data_controlUS=cell2mat(raw_data_controlUS(4:end,:));
    raw_data=cur_nan_col(raw_data);
    raw_data_controlCS=cur_nan_col(raw_data_controlCS);
    raw_data_controlUS=cur_nan_col(raw_data_controlUS);
    %extra data
    A=[];B=[];A_control_CS=[];B_control_CS=[];A_control_US=[];B_control_US=[];
    [A1,B1]=get_X_from_raw(n1,raw_data,rownum_eachfish);
    [A_control_CS1,B_control_CS1]=get_X_from_raw(n1,raw_data_controlCS,rownum_eachfish);
    [A_control_US1,B_control_US1]=get_X_from_raw(n1,raw_data_controlUS,rownum_eachfish);
    [A2,B2]=get_X_from_raw(n2,raw_data,rownum_eachfish);
    [A_control_CS2,B_control_CS2]=get_X_from_raw(n2,raw_data_controlCS,rownum_eachfish);
    [A_control_US2,B_control_US2]=get_X_from_raw(n2,raw_data_controlUS,rownum_eachfish);
    A=get_integrate_X1_X2(A1,A2);B=get_integrate_X1_X2(B1,B2);
    A_control_CS=get_integrate_X1_X2(A_control_CS1,A_control_CS2);B_control_CS=get_integrate_X1_X2(B_control_CS1,B_control_CS2);
    A_control_US=get_integrate_X1_X2(A_control_US1,A_control_US2);B_control_US=get_integrate_X1_X2(B_control_US1,B_control_US2);
    if norm
        B=normlize_zyq_20190320(n,B,normmethod);
        B_control_CS=normlize_zyq_20190320(n,B_control_CS,normmethod);
        B_control_US=normlize_zyq_20190320(n,B_control_US,normmethod);
    end
    if figure_able %画图
        AX=subplot(ceil(max(size(sheet))/3),ceil(max(size(sheet))/ceil(max(size(sheet))/3)),ii);
        Y=[B(:);B_control_CS(:);B_control_US(:)];
        line([2 2],[a b],'linewidth',2,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        line([size(B2,1)+1 size(B2,1)+1],[a b],'linewidth',2,'color',[0.5 0.5 0.5],'linestyle','--');hold on;
        %[min(Y)-0.05*(max(Y)-min(Y)) max(Y)+0.05*(max(Y)-min(Y))]
        %         patch([2 size(B2,1)+1 size(B2,1)+1 2  ],[AX.YLim(1)   AX.YLim(1)  AX.YLim(2) AX.YLim(2)],[0.5 0.5 0.5],...
        %             'edgecolor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.3);hold on;
        if length(n1)==1
            [p1,p2,p3]=plot_mean_area(B,B_control_CS,B_control_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);
            hold on;
            meanB=[];meanB_CS=[];meanB_US=[];
            meanB(:,1)=mean(B(:,1:2:size(B,2)),2,'omitnan');meanB(:,2)=std(B(:,1:2:size(B,2)),[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS(:,1:2:size(B_control_CS,2)),2,'omitnan');meanB_CS(:,2)=std(B_control_CS(:,1:2:size(B_control_CS,2)),[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US(:,1:2:size(B_control_US,2)),2,'omitnan');meanB_US(:,2)=std(B_control_US(:,1:2:size(B_control_US,2)),[],2,'omitnan');
        elseif length(n1)==2 %n为 m and sd
            [p1,p2,p3]=plot_mean_area(B(:,1:2:size(B,2)),B_control_CS(:,1:2:size(B_control_CS,2)),B_control_US(:,1:2:size(B_control_US,2)),...
                sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold on;
            meanB=[];meanB_CS=[];meanB_US=[];
            meanB(:,1)=mean(B(:,1:2:size(B,2)),2,'omitnan');meanB(:,2)=std(B(:,1:2:size(B,2)),[],2,'omitnan');
            meanB_CS(:,1)=mean(B_control_CS(:,1:2:size(B_control_CS,2)),2,'omitnan');meanB_CS(:,2)=std(B_control_CS(:,1:2:size(B_control_CS,2)),[],2,'omitnan');
            meanB_US(:,1)=mean(B_control_US(:,1:2:size(B_control_US,2)),2,'omitnan');meanB_US(:,2)=std(B_control_US(:,1:2:size(B_control_US,2)),[],2,'omitnan');
        end
        for zz=1:length(p1) p1(zz).Color=color_line; end
        if isflash && ~isempty(intersect(n,[1,2,4,5,6,7])) && (size(B,2)==size(name_fishnumber,2) || size(B,2)==size(name_fishnumber,2)*2)
            B(:,1:length(n))=[];
        end
        [p1,p2,p3]=plot_mean_area_errorbar(meanB,meanB_CS,meanB_US,sheet{ii},ii,name_fishnumber,name_fishnumber_CS,name_fishnumber_US);hold on;
        for zz=1:length(p1) p1(zz).Color='r';p1(zz).LineWidth=3;p1(zz).MarkerEdgeColor='r';p1(zz).MarkerFaceColor='r'; end
        text_n(n,B,B_control_CS,B_control_US,figure_able,figure_able_contolCS,figure_able_contolUS);
        AX.XTickLabel =[[]; 'Hab.';string([1:size(B2,1)]');'Tst.'];
        set(gca,'linewidth',2,'fontsize',20);
        ylim([a b]);
    end
end %n为 m and sd

% ind_cut_trial={};
% ind_cut_trial_CS={};
% for ii=1:5
% [behaviorname,behaviorpath]=uigetfile(['E:\A_Data_lightsheet\Data_vmat\20190408' '.mat'],'behavior');
% [fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([behaviorpath,behaviorname]);
% a=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
% ind_cut_trial{ii,1}=unique(a);
% a=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start) & re_startpoint(:,2)<(frameb.us_start-1)),1);
% ind_cut_trial_CS{ii,1}=unique(a);
% end
