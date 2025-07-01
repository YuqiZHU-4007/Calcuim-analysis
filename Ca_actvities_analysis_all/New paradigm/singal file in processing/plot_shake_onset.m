clc;clear all;close all;
[inputname,inputpath]=uigetfile('G:\data\.mat','behavior');
load([inputpath,inputname]);

set(0,'defaultfigurecolor','w');
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([inputpath,inputname]);

h=figure;
title(inputname,'Interpreter','none','fontsize',13);
plot_behavior_onset(delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
saveas(h,[inputpath '\'  'behavior onset' '.tif']);
savefig(h,[inputpath '\'  'behavior onset']);

%Í³¼ÆÍ¼
hh=figure;
plot_behavior_statistics(fs,frame,frameb,trial,re_startpoint);
saveas(hh,[inputpath '\'  'behavior statistics' '.tif']);
savefig(hh,[inputpath '\'  'behavior statistics']);