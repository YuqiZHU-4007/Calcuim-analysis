function yaxis_cut(data)
% ����ضϺ�ͼ��
% ���ߣ���³�¼� - ����԰ http://www.cnblogs.com/kailugaji/
% ����
% ��������
x_min=0; %������̶���Сֵ
x_interval=1; %������̶ȼ������
x_max=30; %������̶����ֵ
y_interval=0.02;  %�����������̶ȼ������
y_max=0.3; %����̶����ֵ
y_break_start=0.05; % �ضϵĿ�ʼֵ
y_break_end=0.15; % �ضϵĽ���ֵ
 
adjust_value=0.4*y_interval; %΢���ضϴ�y����
uptate_num=y_break_end-y_break_start-y_interval; %��ߴ���������ƽ�ƴ�С
 
% �����ضϽ���λ�õ���Щ����ͳͳ����ƽ��uptate_num������
 data(find(data>y_break_end))=data(find(data>y_break_end))-uptate_num;
% �������ߵĸ��������޸ģ�����������4�� 
set(gcf,'color','w') %���汳�����
% ������ض�����
ylimit=get(gca,'ylim');
location_Y=(y_break_start+adjust_value-ylimit(1))/diff(ylimit);
t1=text(0, location_Y,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t1,'rotation',90);
t2=text(1, location_Y,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t2,'rotation',90);
 
% ���¶���������̶�
ytick=0:y_interval:y_max;
set(gca,'ytick',ytick);
ytick(ytick>y_break_start+eps)=ytick(ytick>y_break_start+eps)+uptate_num;
for i=1:length(ytick)
   yticklabel{i}=sprintf('%d',ytick(i));
end
set(gca,'yTickLabel', yticklabel, 'FontSize', 12, 'FontName', 'Times New Roman'); %�޸��������ơ�����
%saveas(gcf,sprintf('Break_Y_Axis.jpg'),'bmp'); %����ͼƬ