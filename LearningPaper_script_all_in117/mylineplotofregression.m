function Formu=mylineplotofregression(x,y,xlab,ylab,tit)
x=x(~isnan(x));y=y(~isnan(y));
p=polyfit(x,y,1);%�?次拟�?;
yfit=polyval(p,x);%求拟合后的y�?;
mdl=fitlm(x, y);
r2 = mdl.Rsquared.Ordinary;%1 - (sum((yfit- y).^2) / sum((y - mean(y)).^2));%即一元线性拟合的R平方
a = num2str(p(1),'%.02f');%即y=ax+b中的a�?
b = num2str(p(2),'%.02f');%即y=ax+b中的b�?
Formu = ['R^2=',num2str(r2,'%.03f')];%这里字符串是拟合公式和R平方
[~,I]=sort(x);
plot(x(I),y(I),'k.',x(I),yfit(I),'k-','Markersize',18,'linewidth',2);%画图;这是散点的大小颜色样�?
%text(max(x),0.45,Formu,'color','r','FontSize',12) %附注说明，即在图上插入一元线性拟合公式和R平方，前两个参数是位置，接着是插入的文本，最后设置字体的大小
%set(gca,'xtick',0.1:0.1:max(x),'XTickLabel',xlab,'xticklabelrotation',45)%这里我设置的是X轴上的数字不可见
%xlim([0 0.8]);ylim([0 0.6]);%ylabel(ylab,'Rotation',90)%这里是y轴的标签
title(tit);box off;
end