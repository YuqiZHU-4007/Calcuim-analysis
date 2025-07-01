function plot_rawtrace_trials(acti,event,fs,frame,trial,startpoint,type)

lab='dfdf';

switch type
case 1%»­hab£¬acq£¬test
    t = tiledlayout('flow');
    set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.65]) ;
    nexttile([1 2]);
    set(gca,'position',[0.01,0.05,0.30,0.18]);
    plot_rawtrace(acti,[trial.spon_bef(3)*frame.per_cycle+1 trial.hab(3)*frame.per_cycle],event,fs,frame,trial,startpoint)
    title(['Pre cond']);ylabel(lab);legend('off');
    nexttile([1 2]);
    set(gca,'position',[0.13,0.13,0.30,0.18]);
    plot_rawtrace(acti,[(trial.acq(3))*frame.per_cycle+1 (trial.test(3))*frame.per_cycle],event,fs,frame,trial,startpoint)
    title(['Post cond']);xlabel('Time(s)');ylabel(lab);legend('off');
    nexttile([1 4]);
    plot_rawtrace(acti,[trial.hab(3)*frame.per_cycle+1 (trial.acq(3))*frame.per_cycle],event,fs,frame,trial,startpoint)
    title(['Cond']);xlabel('Time(s)');ylabel(lab);legend('off');
    
    case 2%»­spon£¨Ç°ºó¶¼»­£©
    set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
    subplot(2,1,1);
    plot_rawtrace(acti,[1 trial.spon_bef(3)*frame.per_cycle],event,fs,frame,trial,startpoint)
    title(['Spon. before habituation']);ylabel(lab);legend('off');
    
    subplot(2,1,2);
    plot_rawtrace(acti,[trial.test(3)*frame.per_cycle+1 trial.total*frame.per_cycle],event,fs,frame,trial,startpoint)
    title(['Spon. after habituation']);ylabel(lab);legend('off');
end
% lab='dfdf';
% set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
% subplot(3,1,1);
% plot_rawtrace(acti(:,ii),[trial.spon_bef(3)*frame.per_cycle+1 trial.hab(3)*frame.per_cycle],event)
% title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Hab']);ylabel(lab);
% subplot(3,1,2);
% plot_rawtrace(acti(:,ii),[trial.hab(3)*frame.per_cycle+1 (trial.acq(3))*frame.per_cycle],event)
% title(['Acq']);xlabel('time(s)');ylabel(lab);
% subplot(3,1,3);
% set(gca,'position',[0.13,0.130780392992875,0.307001594896332,0.18408594969612]);
% plot_rawtrace(acti(:,ii),[(trial.acq(3))*frame.per_cycle+1 (trial.test(3))*frame.per_cycle],event)
% title(['Test']);xlabel('time(s)');ylabel(lab);
% saveas(h2,[outpath '\raw trace\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
% 
% h3=figure;
% aa=acti(:,ii);a=aa([1:trial.spon_bef(3)*frame.per_cycle trial.test(3)*frame.per_cycle+1:trial.total*frame.per_cycle]);
% set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
% subplot(2,1,1),
% plot([1:frame.per_cycle*trial.total]*fs.ca,aa);hold on;xlim([1 trial.spon_bef(3)*frame.per_cycle]*fs.ca);hold on;ylim([min(a(:)) max(a(:))]);
% %plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'g');hold on
% scatter(event.ind*fs.ca,aa(event.ind),9,'k','filled');hold on;
% %behavior
% onset=[];onset=startpoint;
% scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),12,[1 0.5 0],'filled');hold on;
% title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) 'Spon_bef']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
% 
% subplot(2,1,2),
% plot([1:frame.per_cycle*trial.total]*fs.ca,aa);hold on;xlim([trial.test(3)*frame.per_cycle+1 trial.total*frame.per_cycle]*fs.ca);hold on;ylim([min(a(:)) max(a(:))]);
% % plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'g');hold on
% scatter(event.ind*fs.ca,aa(event.ind),9,'k','filled');hold on;
% %behavior
% onset=[];onset=startpoint;
% scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),12,[1 0.5 0],'filled');hold on;
% title(['Spon_aft']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
% saveas(h3,[outpath '\raw trace-sp\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);%ylim([min(a(:))-b max(a(:))])
% close(h3);