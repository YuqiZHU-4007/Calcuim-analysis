colorCS=[0.5 0.5 0.5];
for jj=2
    a=[-0.01 0.03];
    switch jj
        case 1
            x=act_all_hab;isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
        case 2
            x=act_all_acq;isus=true;%set(gcf,'position',[50 50 0.6*k(3) 0.85*k(4)]);
        case 3
            x=trace.test;isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
    end
    for ii=[4 7 10]
        figure,
        ind=find(gIX_iii==ii);
        %step=(1/(length(unique(gIX_iii))+1)*1)/(length(unique(gIX_iii))-1);
        %high=1/(length(unique(gIX_iii))+3);
        %ax=subplot('position',[0.1 0.05+(ii-1)*(1.5*step+high) 0.8 high]);%length(unique(gIX_iii)),1,ii,  
        patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_start:frame.per_cycle:size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        %plot(x(:,cIX_iii(ind)),'linewidth',0.5,'color',[0.5 0.5 0.5]);hold on
        if isus
            l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
        end
        plot(smooth(mean(x(:,cIX_iii(ind)),2)),'linewidth',2,'color',clrmap(ii,:,:));hold on
        ylim(a);box on;
    end
end

