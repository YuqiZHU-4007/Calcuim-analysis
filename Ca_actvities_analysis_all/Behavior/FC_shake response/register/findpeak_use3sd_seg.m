%thr=2*sd
%分段，分成HAB,ACQ,TEST,适用于基线变化的情况
% trial=[0 15 45 60];
% x=delta_r_bef1;
% bin=[10 10 10];
% a=[1 1 1];
function [y,startpoint,thr]=findpeak_use3sd_seg(x,bin,a,trial,totframe_percycle)
%totframe_percycle=1200;
for i=1:length(trial)-1
    %trialnum=trial(i)+1:trial(i+1);
    xx=x((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle); %输入信号
    [yy,~,thr(i)]=findpeak_use3sd(xx,bin(i),a(i),totframe_percycle);
    y((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle)=yy;
end
startpointt=find(y~=0);
startpoint=intersect(startpointt,startpointt(find(y(max(startpointt-1,1))==0)));

figure,
for i=1:length(trial)-1
    plot((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle,x((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle),'b');hold on;
    plot((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle,y((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle),'r');hold on;
    line((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle,thr(i)*ones(size(x((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle))),'color','r');hold on;
    line((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle,-thr(i)*ones(size(x((trial(i))*totframe_percycle+1:trial(i+1)*totframe_percycle))),'color','r');hold on;    
end
scatter(startpoint,y(startpoint),'r');hold off
