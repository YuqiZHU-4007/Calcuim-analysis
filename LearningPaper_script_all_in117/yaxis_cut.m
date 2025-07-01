function yaxis_cut(data)
% 纵轴截断后图像
% 作者：凯鲁嘎吉 - 博客园 http://www.cnblogs.com/kailugaji/
% 数据
% 参数设置
x_min=0; %横坐标刻度最小值
x_interval=1; %横坐标刻度间隔距离
x_max=30; %横坐标刻度最大值
y_interval=0.02;  %纵坐标两个刻度间隔距离
y_max=0.3; %纵轴刻度最大值
y_break_start=0.05; % 截断的开始值
y_break_end=0.15; % 截断的结束值
 
adjust_value=0.4*y_interval; %微调截断处y坐标
uptate_num=y_break_end-y_break_start-y_interval; %最高处曲线向下平移大小
 
% 超过截断结束位置的那些曲线统统向下平移uptate_num个长度
 data(find(data>y_break_end))=data(find(data>y_break_end))-uptate_num;
% 根据曲线的个数进行修改，这里曲线是4条 
set(gcf,'color','w') %后面背景变白
% 纵坐标截断设置
ylimit=get(gca,'ylim');
location_Y=(y_break_start+adjust_value-ylimit(1))/diff(ylimit);
t1=text(0, location_Y,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t1,'rotation',90);
t2=text(1, location_Y,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t2,'rotation',90);
 
% 重新定义纵坐标刻度
ytick=0:y_interval:y_max;
set(gca,'ytick',ytick);
ytick(ytick>y_break_start+eps)=ytick(ytick>y_break_start+eps)+uptate_num;
for i=1:length(ytick)
   yticklabel{i}=sprintf('%d',ytick(i));
end
set(gca,'yTickLabel', yticklabel, 'FontSize', 12, 'FontName', 'Times New Roman'); %修改坐标名称、字体
%saveas(gcf,sprintf('Break_Y_Axis.jpg'),'bmp'); %保存图片