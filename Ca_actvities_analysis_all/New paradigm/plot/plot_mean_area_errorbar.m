
function [p1,p2,p3]=plot_mean_area_errorbar(B,B_control_CS,B_control_US,sheet,num_lab,name_fishnumber,name_fishnumber_CS,name_fishnumber_US)

global norm %logic：是否normalize
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
p1=[];p2=[];p3=[];
if figure_able_contolCS
    p2=errorbar(repmat([1:size(B_control_CS,1)]',1,length(1:2:size(B_control_CS,2))),B_control_CS(:,1:2:size(B_control_CS,2)),...
        B_control_CS(:,1:2:size(B_control_CS,2)),B_control_CS(:,2:2:size(B_control_CS,2)),...
        '-o','color',color_line_controlCS,...
        'LineWidth',linewidth,...
        'MarkerSize',marksize,...
        'MarkerEdgeColor',color_mark_controlCS,...
        'MarkerFaceColor',color_mark_controlCS);hold on;
end
if figure_able_contolUS
    p3=errorbar(repmat([1:size(B_control_US,1)]',1,length(1:2:size(B_control_US,2))),B_control_US(:,1:2:size(B_control_US,2)),...
        B_control_CS(:,1:2:size(B_control_US,2)),B_control_US(:,2:2:size(B_control_US,2)),...
        '-o','color',color_line_controlUS,...
        'LineWidth',linewidth,...
        'MarkerSize',marksize,...
        'MarkerEdgeColor',color_mark_controlUS,...
        'MarkerFaceColor',color_mark_controlUS);hold on;
end
if figure_able_learner
    p1=errorbar(repmat([1:size(B,1)]',1,length(1:2:size(B,2))),B(:,1:2:size(B,2)),...
        B(:,2:2:size(B,2)),...%'color',color_line,...
        '-o',...
        'LineWidth',linewidth,...
        'MarkerSize',marksize,...
        'MarkerEdgeColor',color_mark,...
        'MarkerFaceColor',color_mark);hold on
end
set(gca,'FontSize',20,'xtick',1:size(B,1),'TickLength',[0.01 0.01]);box off;
xlim([0.5 size(B,1)+0.6]);
title(sheet,'FontSize',20);
if norm
    ylabel(['Norm. integrated △F/F0'],'fontsize',15);xlabel('Blocks','fontsize',15);
else
    ylabel(['Integrated △F/F0'],'fontsize',15);xlabel('Blocks','fontsize',15);
end
if lab
    if num_lab==1
        if figure_able_contolCS && figure_able_contolUS && figure_able_learner
            legend([p1 p2 p3],[name_fishnumber,name_fishnumber_CS,name_fishnumber_US]);
        elseif figure_able_contolCS && ~figure_able_contolUS && figure_able_learner
            legend([p1 p2],[name_fishnumber,name_fishnumber_CS]);
        elseif ~figure_able_contolCS && figure_able_contolUS && figure_able_learner
            legend([p1 p3],[name_fishnumber,name_fishnumber_US]);
        elseif  ~figure_able_contolCS && ~figure_able_contolUS && figure_able_learner
            legend([p1],name_fishnumber);
        elseif  ~figure_able_learner
            legend([p2 p3],[name_fishnumber_CS,name_fishnumber_US]);
        end
    end
end
