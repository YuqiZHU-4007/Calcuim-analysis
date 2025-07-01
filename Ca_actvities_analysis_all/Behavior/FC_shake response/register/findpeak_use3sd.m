%thr=2*sd
function [y,startpoint,endpoint,thr]=findpeak_use3sd(x,bin,a,totframe_percycle)
%x=delta_r_bef1; %输入信号
%bin=10; %peak大致宽度，判断该点为peak起始点后，该点到后面bin个点均认为是一个peak，不再判断bin中的点
%a=0.5;%基线大致最大波动,可以略大，但不要大于高峰离群值，基线较小的时候，0.5比较合适
%20181116更新；基线范围由trialnum决定（totframe_percycle改成形参），baseind=[1 min(ceil(trialnum/2),15)]*totframe_percycle;
%20190716更新；加入结束点endpoint,由bin2决定从startpoint后的bin2个点找endpoint
bin2=60;bin3=6;
if nargin<2
    bin=30;
    a=0.5;
end
%totframe_percycle=1200;
trialnum=length(x)/totframe_percycle;
baseind=[1 min(ceil(trialnum/2),15)]*totframe_percycle;%基线的范围，按trial来算
base=x(baseind(1):baseind(2));base(find(abs(base)>a))=mean(base);
%figure,plot(base,'r');hold on;plot(x(5*totframe_percycle+1:15*totframe_percycle));
[pk,loc]=findpeaks(base,'MinPeakHeight',0.1);
thr=mean(abs(pk))+3*std(abs(pk));%mean不行;3*sd有点大，毕竟是起始点（如果要是3*sd，要把峰的宽度减小）
disp(['thr:' num2str(thr)]);
disp(['mean:' num2str(mean(abs(pk)))]);
disp(['sd:' num2str(std(abs(pk)))]);
%x(x<thr & x>0)= thr; x(x>-thr & x<0)=-thr;figure,plot(x,'r');
k=1;y=zeros(size(x));startpointt=[];startpointt(k)=0;
endpointt=[];endpointt(k)=0;
for i=1:length(x)-bin+1
    if i>endpointt(end)
        if abs(x(i))>=thr
            if abs(x(i))>=mean(abs(pk))+4*std(abs(pk)) %适用于变化较大，持续时间短的动
                if sum(abs(x(i+1:i+2))>=thr)>=1
                    startpointt(k)=i;k=k+1;
                    %%find endpoint
                    kk=k-1;
                    t=flip(x(i:min(i+bin2-1,length(x))));endpointt(kk)=startpointt(kk)+bin-1;
                    for jj=1:length(t)-6
                        if abs(t(jj))>=mean(abs(pk))+4*std(abs(pk)) %适用于变化较大，持续时间短的动
                            if sum(abs(t(jj+1:jj+2))>=thr)>=1
                                endpointt(kk)=max(startpointt(kk)+length(t)-jj,startpointt(kk)+bin3-1);
                                y(startpointt(kk):endpointt(kk))=x(startpointt(kk):endpointt(kk));
                                break;
                            end
                        else
                            if sum(abs(t(jj+1:jj+5))>=thr)>=3 %
                                endpointt(kk)=max(startpointt(kk)+length(t)-jj,startpointt(kk)+bin3-1);;
                                y(startpointt(kk):endpointt(kk))=x(startpointt(kk):endpointt(kk));
                                break;
                            end
                        end
                    end %endpoint
                end
            else
                if sum(abs(x(i+1:i+5))>=thr)>=3 %
                    startpointt(k)=i;k=k+1;
                    %%find endpoint
                    kk=k-1;
                    t=flip(x(i:i+bin2-1));endpointt(kk)=startpointt(kk)+bin-1;
                    for jj=1:length(t)-6
                        if abs(t(jj))>=mean(abs(pk))+4*std(abs(pk)) %适用于变化较大，持续时间短的动
                            if sum(abs(t(jj+1:jj+2))>=thr)>=1
                                endpointt(kk)=max(startpointt(kk)+length(t)-jj,startpointt(kk)+bin3-1);
                                y(startpointt(kk):endpointt(kk))=x(startpointt(kk):endpointt(kk));
                                break;
                            end
                        else
                            if sum(abs(t(jj+1:jj+5))>=thr)>=3 %
                                endpointt(kk)=max(startpointt(kk)+length(t)-jj,startpointt(kk)+bin3-1);
                                y(startpointt(kk):endpointt(kk))=x(startpointt(kk):endpointt(kk));
                                break;
                            end
                        end
                    end %endpoint
                end
            end
        end
    end
end
% figure,plot(x);hold on;
% line(1:length(x),thr*ones(size(x)),'color','r');hold on;line(1:length(x),-thr*ones(size(x)),'color','r');hold on;
% plot(y,'r');hold on;
% scatter(startpointt,y(startpointt));hold off
%ind=startpointt(y(startpointt-1)==0);
startpointt=find(y~=0);
startpoint=intersect(startpointt,startpointt(find(y(max(startpointt-1,1))==0)));
endpoint=intersect(startpointt,startpointt(find(y(max(startpointt+1,1))==0)));

%aa=intersect(startpoint,startpoint1);setdiff(startpoint1,aa);
figure,plot(x);hold on;
line(1:length(x),thr*ones(size(x)),'color','r');hold on;line(1:length(x),-thr*ones(size(x)),'color','r');hold on;
plot(y,'r');hold on;
scatter(startpoint,y(startpoint));hold on
scatter(endpoint,y(endpoint),'filled');hold off