
for ii=1:max(size(activities))
    for jj=1:size(a,2)
        a=activities{ii};%_preCS_dfdf
        aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
        %figure,plot_rawtrace_trials(a(:,jj),[],fs,frame,trial,[],1);
        test_smooth(aa',fs,frame,trial)
        test_isoutlier(aa',fs,frame,trial)
    end
end

kk=1;
for ii=1:max(size(activities))
    a=activities{ii};
    for jj=1:size(a,2)
        [~,~,~,error(kk,:)]=test_smooth(rawacti.acq_mean_blocks{ii,jj}(:,3)',fs,frame,trial,polyfit_win,area_win_hab_tst);
        kk=kk+1;
    end
end
figure
boxplot(error,'Notch','on','Labels',{'raw','sgolay','rlowess','movmean'})
title('Compare Random Data from Different Distributions')

for ii=1:4
    subplot(2,2,ii),
    histogram(error(:,ii),50);
    a=sort(error(:,ii));
    Q3=a(ceil(0.75*length(a)));
    Q1=a(ceil(0.25*length(a)));
    QR=Q3-Q1;
    ind1(1:length(find(error(:,ii)>Q3+1.5*QR)),ii)=find(error(:,ii)>Q3+1.5*QR);
    ind2(1:length(find(error(:,ii)<Q1-1.5*QR)),ii)=find(error(:,ii)<Q1-1.5*QR);
    (length(find(error(:,ii)>Q3+1.5*QR))+length(find(error(:,ii)<Q1-1.5*QR)))/size(error,1)
end
kk=0;
for ii=1:max(size(activities))
    a=activities{ii};
    for jj=1:size(a,2)
        kk=kk+1;
        if  ~isempty(find(kk==ind1(:)))
            [~,~,~,error(kk,:)]=test_smooth(rawacti.hab{ii,1}(:,jj)',fs,frame,trial,polyfit_win,area_win_hab_tst);
        end
    end
end

layer=13;num=20;
ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
linewidth=3;legend_location='northwest';
color=[0.5 0.5 0.5];%CS bar color
ColorOrder=[0,0,0;...
    0.3,0.75,0.93;...
    0,0.45,0.74;...
    0,0,1;...
    0.93,0.69,0.13;...
    0.85,0.33,0.10;...
    1,0,0;...
    0.64,0.08,0.18;...
    0.49,0.18,0.56;...
    1,0,1;];
activities={};
for ii = 1:size(activities_new,1)
    if ~isempty(activities_new{ii,1}(:))
        activities{ii,1}=[activities_new{ii,1}(1,:);activities_new{ii,1};activities_new{ii,1}(size(activities_new{ii,1},1),:);];
    else
        activities{ii,1}=[];
    end
end
figure,
plot(activities{layer,1}(:,num))
%plot_rawtrace_trials(activities{layer,1}(:,num),[],fs,frame,trial,startpoint,1);
%area_mean_across_rois=plot_meantrace_acq(activities_preCS_dfdf_aftcorrect{layer,1}(:,num),fs,trial,frame,false,[]);

figure,
for ii=1:2
    subplot(1,2,ii, 'ColorOrder',ColorOrder),
    edge=[300 -100];
    %CS bar
    patch1=patch(([frame.cs_start...
        frame.cs_end...
        frame.cs_end'...
        frame.cs_start]-frame.cs_start)*fs.ca,...
        [edge(1) edge(1) edge(2) edge(2)],...
        color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
    %US line
    l1=line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[edge(1) edge(2)],'color',[1 0 0],'linestyle','--','linewidth',linewidth);hold on
    switch ii
        case 1
            acti=activities_preCS_dfdf_aftcorrect{layer,1}(:,num);
            for ii=1:trial.acq_block_num
                ind1=trial.acq_block_trial*(ii-1)+1+trial.hab(3);
                ind2=trial.acq_block_trial*(ii-1)+1+trial.hab(3)+trial.acq_block_trial-1;
                %subplot(trial.acq_block_num,1,ii),
                %plot(aa((ind1-1)*frame.per_cycle+1:ind1*frame.per_cycle));hold on
                a=mean(reshape(acti((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)',1);
                edge(1)=min(edge(1),min(a));edge(2)=max(edge(2),max(a));
                plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,a,...
                    'linewidth',linewidth);hold on
                legend(['CS';'US';num2str([1:trial.acq_block_num]','%02d')],'location',legend_location);
                title(['Mean across trials of each block:' num2str(layer) '-'  num2str(num)],'FontSize',20,'interprete','none');
                xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
                ylabel('△F/F0');xlabel('Time(s)');
                set(gca,'FontSize',20);box off;
            end
            
        case 2
            acti=rawacti.acq_mean_blocks{layer,num};
            plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
                'linewidth',linewidth);hold on
            edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
            %             acti=rawacti.acq_1st_block{layer,num};
            %             plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
            %                 'linewidth',linewidth);hold on
            %             edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
            legend(['CS';'US';num2str([1:trial.acq_block_num]','%02d')],'location',legend_location);
            title(['Mean across trials of each block:' num2str(layer) '-'  num2str(num)],'FontSize',20,'interprete','none');
            xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
            ylabel('△F/F0');xlabel('Time(s)');
            set(gca,'FontSize',20);box off;
    end
    % edge=[edge(1) edge(2)]+[-0.25 0.25]*(edge(2)-edge(1));
end
figure,set(gca,'ColorOrder',ColorOrder);
for ii=1
    edge=[100 -100];
    %CS bar
    patch1=patch(([frame.cs_start...
        frame.cs_end...
        frame.cs_end'...
        frame.cs_start]-frame.cs_start)*fs.ca,...
        [edge(1) edge(1) edge(2) edge(2)],...
        color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
    %US line
    l1=line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[edge(1) edge(2)],'color',[1 0 0],'linestyle','--','linewidth',linewidth);hold on
    acti=rawacti.acq_mean_blocks{layer,num};
    plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
        'linewidth',linewidth);hold on
    edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
    %             acti=rawacti.acq_1st_block{layer,num};
    %             plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
    %                 'linewidth',linewidth);hold on
    %             edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
    legend(['CS';'US';num2str([1:trial.acq_block_num]','%02d')],'location',legend_location);
    title(['Mean across trials of each block:' num2str(layer) '-'  num2str(num)],'FontSize',20,'interprete','none');
    xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
    ylabel('△F/F0');xlabel('Time(s)');
    set(gca,'FontSize',20);box off;
end

area_win_acq=frame.cs_start:frame.us_start-1;
area_win_hab_tst=frame.cs_start:frame.cs_end-1;
ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;%取拟合时间窗
type={'lowest','lowest_all','first','mean','sum','event'};
figure,
for jj=1:size(type,2)
    area=[];
    area=calculate_integrate_dfdf(rawacti.acq_mean_blocks{layer,num},area_win_acq,type{jj},ref_win);
    subplot(3,2,jj)
    x = 1:trial.acq_block_num;
    plot(x,area,'-o',...
        'LineWidth',3,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
    ylabel('Integrated △F/F0');xlabel('Blocks');
    title(type{jj},'FontSize',20,'interprete','none');
    %ylim([-0.1 0.05])
end

figure,set(gca,'ColorOrder',ColorOrder);
for ii=1
    edge=[100 -100];
    %CS bar
    patch1=patch(([frame.cs_start...
        frame.cs_end...
        frame.cs_end'...
        frame.cs_start]-frame.cs_start)*fs.ca,...
        [edge(1) edge(1) edge(2) edge(2)],...
        color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
    %US line
    l1=line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[edge(1) edge(2)],'color',[1 0 0],'linestyle','--','linewidth',linewidth);hold on
    acti=rawacti.hab{layer}(:,num);
    plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
        'linewidth',linewidth);hold on
    acti=rawacti.tst{layer}(:,num);
    plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
        'linewidth',linewidth);hold on
    edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
    %             acti=rawacti.acq_1st_block{layer,num};
    %             plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,acti,...
    %                 'linewidth',linewidth);hold on
    %             edge(1)=min(edge(1),min(acti(:)));edge(2)=max(edge(2),max(acti(:)));
    %legend(['CS';'US';num2str([1:trial.acq_block_num]','%02d')],'location',legend_location);
   % title(['Mean across trials of each block:' num2str(layer) '-'  num2str(num)],'FontSize',20,'interprete','none');
    xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
    ylabel('△F/F0');xlabel('Time(s)');
    set(gca,'FontSize',20);box off;
end
figure,
for jj=1:size(type,2)
    area=[];
    area=calculate_integrate_dfdf([rawacti.hab{layer}(:,num),rawacti.tst{layer}(:,num)],area_win_hab_tst,type{jj},ref_win);
    %area=calculate_integrate_dfdf(rawacti.tst{layer}(:,num),area_win_hab_tst,type{jj},ref_win);
    subplot(3,2,jj)
    x = 1:2;
    plot(x,area,'-o',...
        'LineWidth',3,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
    ylabel('Integrated △F/F0');xlabel('Blocks');
    title(type{jj},'FontSize',20,'interprete','none');
    %ylim([-0.1 0.05])
end
% figure,
% plot(x,area.CS_acq_block{layer,1}(:,num),'-o',...
%         'LineWidth',3,...
%         'MarkerEdgeColor','r',...
%         'MarkerFaceColor','r')
%     set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
%     ylabel('Integrated △F/F0');xlabel('Blocks');
%     title(type{jj},'FontSize',20);