function plot_test_1(x,B,cIX2,gIX2,clrmap,a,frame,fs,type)

colnum=5;
switch type
    case 1 %画在一起
        for i = 1:length(unique(gIX2))
            subplot(ceil(length(unique(gIX2))/colnum),colnum,i)
            clust = find(gIX2==i);
            %     plot(x,B(:,cIX2(clust)),'color',[0.5 0.5 0.5],...
            %         'LineWidth',1,...
            %         'markersize',3);
            %     hold on
            plot(x,mean(B(:,cIX2(clust)),2),'-o','color',clrmap(i,:),...
                'LineWidth',2,...
                'markersize',3);
            hold on
            if ~exist('a') || isempty(a)
                a=[min(mean(B(:,cIX2(clust)),2)) max(mean(B(:,cIX2(clust)),2))];
            end
            %     line(([frame.cs_start frame.cs_start]-frame.cs_start)*fs,a,'color','b','linewidth',1.5,'linestyle','--');hold on;
            %     line(([frame.cs_end frame.cs_end]-frame.cs_start)*fs,a,'color','b','linewidth',1.5,'linestyle','--');hold on;
            set(gca,'FontSize',12);box off;%ylim([0 1.5]);
            %set(gca,'visible','off');
            %text(5,(min(a)+max(a))/2,num2str(length(clust)),'fontsize',25);hold on;
            %ylabel('△F/F0');
            %xlabel('Time');
            box on;
            ylim(a);
            xlim([0.5 length(x)+0.5]);
        end
        hold off
    case 2 %前后一个点分开
        for i = 1:length(unique(gIX2))
            AX=subplot(ceil(length(unique(gIX2))/colnum),colnum,i);
            clust = find(gIX2==i);
            a=mean(B(:,cIX2(clust)),2);a=(a-a(1))/a(1);
            plot(x(2:end-1),a(2:end-1)-(a(2)),'-o','color',clrmap(i,:),...
                'LineWidth',2,...
                'markersize',3);
            hold on
            scatter([x(1) x(end)],[a(1) a(end)],30,clrmap(i,:),'filled');hold on;
            if ~exist('a') || isempty(a)
                a=[min(mean(B(:,cIX2(clust)),2)) max(mean(B(:,cIX2(clust)),2))];
            end
            set(gca,'FontSize',12);box off;%ylim([0 1.5]);
            box on;
            %ylim(a);
            xticks(x)
            AX.XTickLabel =[[],'Hab.';string(x(2:end-1)'-1);'Tst.'];
            xlim([0.5 length(x)+0.5]);set(gca,'xcolor','w','ycolor','w');
        end
        
        hold off
end
