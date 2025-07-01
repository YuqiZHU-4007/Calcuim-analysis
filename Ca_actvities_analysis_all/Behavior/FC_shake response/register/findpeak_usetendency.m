%三角形模拟
%20181121改：加bin参数，起始点间距小于bin算做一次
function [y,startpoint,thr]=findpeak_usetendency(x1,thr,bin)
x=abs(x1);%figure,plot(x)
if nargin<2
%thr=0.2;
[pk,loc]=findpeaks(x);
pk=sort(pk);thr=pk(fix(length(pk)*0.95));disp(['thr:' num2str(thr)]);
bin=10;
end
xx= smoothdata(x,'movmean',3);%这个效果较好,但是边界不太好
%xx=x;
%xx(x>=thr & [x(1) x(1:end-1)]<thr & [x(2:end) x(end)]<thr)=0;
%figure,plot(xx,'r');hold on;plot(x);
%xx(x<thr)=0;figure,plot(xx)
y=zeros(size(x));
for i=2:length(xx)-2
    if xx(i)>=thr 
        if xx(i-1)<thr && xx(i+1)<thr && xx(i+2)<thr
        y(i)=0;
        elseif xx(i-1)<thr && xx(i+1)>=thr && xx(i+2)>=thr
            y(i)=x(i);
        elseif xx(i-1)>thr &&  (xx(i-2)>thr || xx(i+1)>=thr)
            y(i)=x(i);
        end
    elseif xx(i)<thr && xx(i-1)>=thr && xx(i+1)>=thr %动的过程中可能会有较低值
        y(i)=x(i);
    else
        y(i)=0;
    end
end
startpointt1=find(y~=0);
startpointt=intersect(startpointt1,startpointt1(find(y(max(startpointt1-1,1))==0)));
%startpoint(find(startpoint(1:end-1)-startpoint(2:end)>-bin)+1)=[];%20181121改：加bin参数，起始点间距小于bin算做一次
startpoint=startpointt;
for ii=1:length(startpointt)-1
    if ii<=length(startpoint)-1
        kk=ii;
        while startpoint(kk+1)<=startpoint(kk)+bin-1
        %startpoint(ii)-startpoint(ii+1)>-bin
            startpoint(kk+1)=[];kk=min(kk,length(startpoint)-1);
        end
    end
end
%startpoint(find(startpoint==0))=[];

figure,plot(x);hold on;
line(1:length(x),thr*ones(size(x)),'color','r');hold on;line(1:length(x),-thr*ones(size(x)),'color','r');hold on;
plot(y,'r');hold on;
scatter(startpoint,y(startpoint));hold off