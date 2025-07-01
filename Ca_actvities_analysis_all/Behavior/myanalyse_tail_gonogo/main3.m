clear all;
close all;
[filename, pathname] = uigetfile( {'*.mat'}, 'Pick a file');
% filename='alltheta.mat';
% pathname='C:\light data\20170929\fish1\';
% Taildata=load(filename,pathname);
load([pathname filename]);seg='10.9,fish2';
% filepath =[dirpath filename];
%savepath = checkpath([dirpath '/' imgname '_find_tail/']);
savepath = pathname;
set(0, 'DefaultFigureColor', 'white');

% theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
% theta(theta > 180) = theta(theta > 180) - 360;
% theta(theta < -180) = theta(theta < -180) + 360;
% theta = -theta;
% alltheta = 0 * ones(nframes, 1);
% alltheta(movind) = theta;
%%%%%%
fs=240;
fs_ca=25;
Interval = 50;Latency = 30;
trailnum=floor(length(alltheta)/12000)+1;
sti_start=7105;duration=240*1;
tindex=5;bin1=81;bin2=80;
%tAvg = (1/fs-Latency:1/fs:Interval-Latency)';
time=1/fs:1/fs:Interval;
time_sti=30:1/fs:31-1/fs;
time_befsti=1:1/fs:Latency-1/fs;
%%%%找所有动 和刺激画在一起
if length(alltheta)<trailnum*12000
alltheta(length(alltheta)+1:12000*trailnum)=0;
end
%%%去前后
[taildeg_everytrial,taildeg_everytrial1]=mov_detection_1(alltheta,bin1,bin2,trailnum);
%%%去小尾动%%%
t_dur=zeros;
t_dur=find_t_dur(taildeg_everytrial1,trailnum);
base1=20;%%高峰与干扰的均值差别
base2=0.45*10;%%%去峰值前端 dif的值 春峰师兄的观察值
base3=30;%%%去双峰之间的 detection距离
[t_dur_cell,taildeg_everytrial1]=move_detection_3(t_dur,taildeg_everytrial1,trailnum,base1,base2,base3,bin2);%%如果有很长的连续动，加限定条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:trailnum
%     a=t_dur(:,1,i);%a(find(a==0),:)=[];
%     b=t_dur(:,2,i);%b(find(b==0),:)=[];
%     for j=1:length(a)
%         ind=zeros;
%         if a(j)~=0
%             text=taildeg_everytrial1(i,a(j):a(j)+b(j)-1).*10;
%             dif=abs(diff(text));
%             ind=find(abs(diff(text))>=50);%%50 认为是高峰 与低尾动干扰的阈值
%             if (sum(abs(text))/length(text)<=base1 && length(ind)==0 ) %%均值小于阈值，且持续时间长于120&& b(j)>=120 判断正负干扰全在一侧
%                 t_dur(j,:,i)=0;taildeg_everytrial1(i,a(j):a(j)+b(j)-1)=0;%%去干扰         
%             elseif (dif(1)<=base2  && length(ind)~=0) %%%去峰值前端   
%                 ind=find([dif(2:end) dif(end)]>=base2 &  dif>=base2 & [dif(1) dif(1:end-1)]<base2);%[dif(3:end) dif(end) dif(end)]>=base2 &
%                 %taildeg_everytrial1(i,a(j):a(j)+ind(1))=0;
%                 if length(ind)~=0
%                 ind(1)=ind(1)+a(j);
%                 t_dur(j,:,i)=[ind(1) b(j)-(ind(1)-a(j)+1)];
%                 end
%             end
%             if (b(j)>200 && length(ind)~=0) %%%去双峰之间
%                 ind=find(abs(diff(text))>=45);
%                 mind=find([ind(2:end) ind(end)]-ind>=base3);%%%峰值连续变化index的插值的距离，大于30 认为有第二次动
%                 l=length(mind);
%                 if l~=0
%                 if l>1
%                     for z=1:l
%                     mind(z)=ind(mind(z)+1)+a(j)-1;
%                     taildeg_everytrial1(i,mind(z)-bin2:mind(z)-1)=0;
%                     end
%                     t_dur(j,2,i)=mind(1)-bin2-t_dur(j,1,i)+1;
%                     t_dur(size(t_dur,1):size(t_dur,1)+l-1,:,i)=[mind' (a(j)+b(j)-mind)'];
%                 else
%                     mind=ind(mind+1)+a(j)-1;
%                     taildeg_everytrial1(i,mind-bin2:mind-1)=0;
%                     t_dur(j,2,i)=mind-bin2-t_dur(j,1,i)+1;
%                     t_dur(size(t_dur,1)+1,:,i)=[mind a(j)+b(j)-mind];
%                 end
%                 end
%             end
%         end
%     end
%     tt=sortrows(t_dur(:,:,i),1);
%     tt(find(tt(:,1)==0),:)=[];t_dur_cell{i,1}=tt;    
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t_dur(find(t_dur(:,1,:)==0),:,:)=[];
% win=1:8;
% base1=3;base2=2.5;
% [taildeg_everytrial2]=mov_detection_2(taildeg_everytrial1,win,base1,base2,trailnum);
i=21;
% diff(taildeg_everytrial);
% uu=unique(taildeg_everytrial(i,3050:3150));
figure,%ylim([-5,5]);
title(['trial number 5 in ' seg])
subplot(2,1,1),plot((taildeg_everytrial(i,:)));title(['original data-trial number' i  'in ' seg]);
subplot(2,1,2),plot(1:12000,taildeg_everytrial1(i,:),'r');title('process data');%,taildeg_everytrial2(i,:)
saveas(gcf,[savepath,'1','.png']);
% sum(abs(taildeg_everytrial(i,1782:1833)).*100)/51
% sum(abs(taildeg_everytrial(i,3799:3863)).*100)/64
% sum(abs(taildeg_everytrial(i,6946:7235)).*100)/289
% sum(abs(taildeg_everytrial(i,1117:1217)).*10)/100
% t_dur=zeros;
% t_dur=find_t_dur(taildeg_everytrial2,trailnum);
% ttotal_dur=find_t_durn(taildeg_everytrial2,trailnum);
% 
%subplot(3,1,3),plot(taildeg_everytrial2(i,:));
% hfig=0:0.05:trailnum*0.05;
% figure,
% for i=1:trailnum
% subplot(i,1,1);
% plot(time,taildeg_everytrial(i,:));hold on;
% set(gca,'position',[0.1 hfig(i) 0.5 0.5]);
% set(gca,'visible','off');
% % end
% i=1;N=12000;
% win=zeros(1,N);win(1)=1;win(3*fs:1:8*fs)=1;
% a=abs(fft(taildeg_everytrial(1,:)));a=a(1:6000);
% b=abs(fft(taildeg_everytrial(2,:)));b=b(1:6000);
% c=abs(fft(taildeg_everytrial(3,:)));c=c(1:6000);
% t=(0:(6000-1))/fs;
% figure,%line([1 1],[0 1],'color','r');hold on;line([8 8],[0 1],'color','r');hold on;
% subplot(3,1,1);plot(t,a*2/N);title(['fft trial number 1 in ' seg]);line([3 3],[0 1],'color','r');hold on;line([8 8],[0 1],'color','r');hold on;
% subplot(3,1,2);plot(t,b*2/N);title(['fft trial number 2 in ' seg]);line([3 3],[0 1],'color','r');hold on;line([8 8],[0 1],'color','r');hold on;
% subplot(3,1,3);plot(t,c*2/N);title(['fft trial number 3 in ' seg]);line([3 3],[0 1],'color','r');hold on;line([8 8],[0 1],'color','r');hold on;
% t=1:12000;
% figure,subplot(3,1,1);plot(t,real(ifft((fft(taildeg_everytrial(1,:))).*win)));title(['ifft trial number 1 in ' seg]);
% subplot(3,1,2);plot(t,real(ifft((fft(taildeg_everytrial(2,:))).*win)));title(['ifft trial number 2 in ' seg]);
% subplot(3,1,3);plot(t,real(ifft((fft(taildeg_everytrial(3,:))).*win)));title(['ifft trial number 3 in ' seg]);
% 
% figure,
% subplot(2,1,1);plot(t,taildeg_everytrial(1,:));title(['rewdata trial number 1 in ' seg]);
% subplot(2,1,2);plot(t,real(ifft((fft(taildeg_everytrial(1,:))).*win)));title(['ifft trial number 1 in ' seg]);
% figure,
% subplot(2,1,1);plot(t,taildeg_everytrial(2,:));title(['rewdata trial number 1 in ' seg]);
% subplot(2,1,2);plot(t,real(ifft((fft(taildeg_everytrial(2,:))).*win)));title(['ifft trial number 2 in ' seg]);

h=figure;
title(seg);
interval=4;
line([sti_start sti_start],[0 interval*(trailnum+1) ],'Color','r');hold on
line([sti_start+duration sti_start+duration],[0 interval*(trailnum+1) ],'Color','r');hold on
for i=1:42%trailnum
    a=t_dur_cell{i,1}; 
    y=ones(length(a(:,1)),1).*interval*(i-1);
    scatter(a(:,1),y,10,'filled','r');hold on
    scatter(a(:,3),y,10,'filled','k');hold on
    %line([1 12000],[interval*(i-1) interval*(i-1)]);hold on 
    plot(1:12000,taildeg_everytrial(i,:)/10+interval*(i-1));hold on
    xlabel('frame');
    ylabel('tail point 10 degree for every trial');
end
hold off;
saveas(gcf,[savepath,'2','.png']);
 %print(h,'-dpng',[savepath 'total.jpg']); 
 trail_mov_number=zeros;trail_nomov_number=zeros;trail_false_number=zeros;
[trail_mov_number,trail_nomov_number,trail_false_number,sti_tdur,base_tdur]=find_mov_trailnum(t_dur_cell,sti_start,duration,trailnum);
save([savepath '/t-dur.mat'], 'taildeg_everytrial','t_dur_cell', 'sti_tdur','base_tdur','trail_mov_number','trail_nomov_number','trail_false_number','sti_start','duration');
%%%%%%%%%%统计图%%%%%%%%%%%
if trail_false_number==0
    A = [zeros(length(trail_nomov_number),1)' ones(length(trail_mov_number),1)'];
    C = categorical(A,[1 0],{'mov','nonmov'});
else
A = [ones(length(trail_mov_number),1)' zeros(length(trail_nomov_number)+length(trail_false_number),1)'];%NaN(length(trail_false_number),1)'
C = categorical(A,[1 0],{'mov','nomov'}); %NaN ,'false'
end

figure,h=histogram(C,'BarWidth',0.5);title(['mov and nonmov trial number for stimulate' seg]);
saveas(gcf,[savepath,'3','.png']);
% n=[length(trail_nomov_number) length(trail_mov_number) length(trail_false_number)];
% str= mat2cell(num2str(n'));
% text(h,n+0.3,str,'Color','k');
%%%刺激前后动的均值以及持续时间均值，平均动的次数
n=1000;a=zeros(1,n);b=0;tbaseave=0;tstiave=0;
for i=1:trailnum
    z=base_tdur(:,1,i);z(find(z==0),:)=[];
    for j=1:size(base_tdur(:,1,i),1)
    a=a+[taildeg_everytrial(i,(base_tdur(j,1,i):base_tdur(j,1,i)+base_tdur(j,2,i)-1)) zeros(n-base_tdur(j,2,i),1)'];
    tbaseave=tbaseave+base_tdur(j,2,i);
    end
    b=b+size(z,1);
end
c=zeros(1,n);d=0;
for i=1:trailnum
    z=sti_tdur(:,1,i);z(find(z==0),:)=[];
    for j=1:size(sti_tdur(:,1,i),1)
    c=c+[taildeg_everytrial(i,(sti_tdur(j,1,i):sti_tdur(j,1,i)+sti_tdur(j,2,i)-1)) zeros(n-sti_tdur(j,2,i),1)'];
    tstiave=tstiave+sti_tdur(j,2,i);
    end
    d=d+size(z,1);
end
figure,subplot(2,1,1);plot(a/b);title(['average all baseline tail movement' seg]);
subplot(2,1,2);plot(c/d);title('average all stimulating tail movement');
saveas(gcf,[savepath,'4','.png']);
tt=zeros;
for i=1:trailnum
    t=t_dur_cell{i,1};l=length(tt);
    tt(l+1:size(t,1)+l)=t(:,1);
end
figure,hist(tt,100);title(['all movement' seg]);hold on
line([sti_start sti_start],[0 30 ],'Color','r');hold on
line([sti_start+duration sti_start+duration],[0 30 ],'Color','r');
tbaseave/b
tstiave/d
saveas(gcf,[savepath,'5','.png']);
% A = [zeros(tbaseave/b,1)' ones(tstiave/d,1)'];
% C = categorical(A,[1 0],{'baseline movement','stimulating movement'});
% figure,h=histogram(C,'BarWidth',0.5);title('durarion average-9.29,fish1');