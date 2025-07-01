%detect shake 2018/08/23
figure,plot(diff(delta_r_bef1));hold on
plot(delta_r_bef1,'r')

xx=x;xx(abs(xx)<mean(x))=0;y=zeros(size(x));
for i=1:length(xx)-2
    if sum(abs(xx(i:i+2))>thr)>=3
        y(i:i+2)=x(i:i+2);
    end
end
figure,plot(x);hold on;
line(1:length(x),thr*ones(size(x)),'color','r');hold on;line(1:length(x),-thr*ones(size(x)),'color','r');hold on;
plot(y,'r');hold on;
%scatter(startpoint,y(startpoint));hold off
%峰值
p=0.98;
x=delta_r_bef1;
[pk,loc]=findpeaks(x(1:1200));pk=sort(pk);%1.37%为真实峰值
thr=max(pk(1:fix(length(pk)*p)));%去掉90%的值
disp(num2str(thr));disp(length(find(pk>thr)));
x(find(abs(x)<thr))=0;
y=zeros(size(x));k=1;bin=10;
for jj=1:size(x,2)-bin
    if sum(y(jj:jj+bin-1)==zeros(size(y(jj:jj+bin-1))))==10 %startpoint(end)+bin-1<jj
        if abs(x(jj))>=thr && sum([abs(x(jj+1:jj+3))>=thr])>=2
            y(jj:jj+bin-1)=x(jj:jj+bin-1);startpointt(k)=jj;
            k=k+1;
        end
    end  
end
figure,plot(delta_r_bef1);hold on
plot(y,'r');hold on;

%各种滤波（中值，小波，高斯）
x=delta_r_bef1;
y=medfilt1(x,30);figure,plot(delta_r_bef1);hold on;plot(y,'r');hold on;


[pk,loc]=findpeaks(x);
thr=mean(pk);
x(find(abs(x)<=thr))=0;
figure,plot(delta_r_bef1);hold on
plot(y,'r');hold on;
scatter(loc,pk);