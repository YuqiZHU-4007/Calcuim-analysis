function plot_shakenum(re_startpoint,trials,frameb,fs_behavior,bin)
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trials(2) & re_startpoint(:,1)<=trials(3)),2);xx=[];
xx(1)=length(find(x<frameb.cs_start & x>=frameb.cs_start-4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.cs_start-1);
xx(2)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs_behavior;%(frameb.cs_end-frameb.cs_start);
xx(3)=length(find(x>=frameb.cs_end & x<frameb.cs_end+4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.per_cycle-frameb.cs_end)+1;
summ=trials(1)*bin;
histogram('Categories',{'Bef CS(4.8s)','CS','Aft CS(4.8s)'},'BinCounts',[xx(1),xx(2),xx(3)]/summ);
higher_than_edge=max(xx)/summ*0.03;
%Before CS
f= table([xx(2);trials(1)*bin-xx(2)],[xx(1);trials(1)*bin-xx(1)],'VariableNames',{'CS','NoCS'},'RowNames',{'MOV','NoMOV'});f;
[h,p,stats]=plot_fishertest_p(f,xx(1)/summ,xx(2)/summ,categorical({'Bef CS(4.8s)','CS'}),higher_than_edge);
%After CS
%higher_than_edge=higher_than_edge*2;
f= table([xx(2);trials(1)*bin-xx(2)],[xx(3);trials(1)*bin-xx(3)],'VariableNames',{'CS','NoCS'},'RowNames',{'MOV','NoMOV'});f;
if max(xx(2),xx(3))==max(xx(1),xx(2))
    higher_than_edge=higher_than_edge*2;
end
[h,p,stats]=plot_fishertest_p(f,xx(2)/summ,xx(3)/summ,categorical({'CS','Aft CS(4.8s)'}),higher_than_edge);
ylim([0 max(max(xx),0.0001)/summ+higher_than_edge*3]);
set(gca,'FontSize',13);
% onset=[];onset=re_startpoint;%行为
% x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2) & re_startpoint(:,1)<=trial.hab(3)),2);xx=[];
% xx(1)=length(find(x<frameb.cs_start & x>=frameb.cs_start-4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.cs_start-1);
% xx(2)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs_behavior;%(frameb.cs_end-frameb.cs_start);
% xx(3)=length(find(x>=frameb.cs_end & x<frameb.cs_end+4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.per_cycle-frameb.cs_end)+1;
% histogram('Categories',{'Bef CS(4.8s)','CS','Aft CS(4.8s)'},'BinCounts',[xx(1),xx(2),xx(3)]);
% %Before CS
% f= table([xx(2);trial.hab(1)*bin2-xx(2)],[xx(1);trial.hab(1)*bin2-xx(1)],'VariableNames',{'CS','NoCS'},'RowNames',{'MOV','NoMOV'});f;
% [h,p,stats]=plot_fishertest_p(f,xx(1),xx(2),categorical({'Bef CS(4.8s)','CS'}),higher_than_edge);
% %After CS
% %higher_than_edge=higher_than_edge*2;
% f= table([xx(2);trial.hab(1)*bin2-xx(2)],[xx(3);trial.hab(1)*bin2-xx(3)],'VariableNames',{'CS','NoCS'},'RowNames',{'MOV','NoMOV'});f;
% [h,p,stats]=plot_fishertest_p(f,xx(2),xx(3),categorical({'CS','Aft CS(4.8s)'}),higher_than_edge);
% ylim([0 max([xx(3),xx(2)])+higher_than_edge*3]);

% x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2) & re_startpoint(:,1)<=trial.hab(3)),2);xx=[];
% xx{1,1}=x(find(x<frameb.cs_start & x>=frameb.cs_start-4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.cs_start-1);
% h=histogram(xx{1,1},bin,'facecolor','k','Normalization',method,'Visible','off' );
% xx{1,2}(1)=length(find(h.BinCounts==1));xx{1,2}(2)=length(find(h.BinCounts==0));
% 
% xx{2,1}=x(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs_behavior;%(frameb.cs_end-frameb.cs_start);
% h=histogram(xx{2,1},bin,'facecolor','k','Normalization',method,'Visible','off' );
% xx{2,2}(1)=length(find(h.BinCounts==1));xx{2,2}(2)=length(find(h.BinCounts==0));
% 
% xx{3,1}=x(find(x>=frameb.cs_end & x<frameb.cs_end+4.8*fs_behavior));%/4.8*fs_behavior;%(frameb.per_cycle-frameb.cs_end)+1;
% h=histogram(xx{3,1},bin,'facecolor','k','Normalization',method,'Visible','off' );
% xx{3,2}(1)=length(find(h.BinCounts==1));xx{3,2}(2)=length(find(h.BinCounts==0));