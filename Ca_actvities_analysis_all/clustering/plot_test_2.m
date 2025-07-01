function [h,ratio,trace_for_ttest]=plot_test_2(A,cIX,gIX,frame,clrmap,type,yli,linevisible,isplotback)
if  ~exist('type','var')
    type=1;
end
if length(size(A))<3
    A=reshape(A,[size(A,1),1,size(A,2)]);
end
trace_for_ttest={};ratio=[];
ColorOrder=[0,0,0;...
    0,1,0;...
    1,0,0;...
    0,0,1;...
    1,1,0;];
set(groot,'defaultAxesColorOrder',ColorOrder)
h=figure('position',[0,0,1910,1000]);
%m=[-0.04 0.12];
switch type
    case 1
        if linevisible
            l3=plot(repmat([(frame.us_start):frame.per_cycle:frame.per_cycle*size(A,2)],2,1)',[min(yli(:)) max(yli(:))],'r','LineStyle','-','linewidth',1.2);hold on
        end
        %CS bar
        color=[0.5 0.5 0.5];
        patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]'...
            [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
            [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
            [frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]']',...
            repmat([min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],size(A,2),1)',...
            color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
        ratio=[];
        if isplotback
            plot(mean(reshape(A(:,:,:),frame.per_cycle*size(A,2),[]),2),'linewidth',1,'color',[0.5 0.5 0.5]);hold on
        end
        for kk=unique(gIX)'
            ind=cIX(find(gIX==kk)); ratio(kk,1)=length(ind);
            m=mean(reshape(A(:,:,ind),frame.per_cycle*size(A,2),[]),2);
            plot(m,'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1); hold on
            shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1); hold on
            %plot(reshape(A(:,:,ind),frame.per_cycle*size(A,2),[]),'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            %plot_rawtrace(m,[1:length(m)],[],fs,frame,trial,[]);hold on;
        end
        ylim(yli);
    case 2
        colnum=3;ratio=[];
        for kk=unique(gIX)'
            subplot(ceil(max(unique(gIX))/colnum),colnum,kk);
            if linevisible
                l3=plot(repmat([(frame.us_start):frame.per_cycle:frame.per_cycle*size(A,2)],2,1)',[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',1.2);hold on
            end
            %CS bar
            color=[0.5 0.5 0.5];
            patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]'...
                [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
                [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
                [frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]']',...
                repmat([min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],size(A,2),1)',...
                color,'edgecolor',color,'facecolor',color,'edgealpha',0.4,'facealpha',0.45);hold on
            if isplotback
                plot(mean(reshape(A(:,:,:),frame.per_cycle*size(A,2),[]),2),'LineStyle','--','linewidth',1.5,'color',[0.5 0.5 0.5]);hold on
            end
            ind=cIX(find(gIX==kk)); ratio(kk,1)=length(ind);
            m=mean(reshape(A(:,:,ind),frame.per_cycle*size(A,2),[]),2);
            plot(m,'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            trace_for_ttest{kk}=reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])';
            %yli=[0 0.1];
            ylim(yli);title([num2str(kk,'%02d')]);
            set(gca,'xcolor','w','ycolor','w');
        end
    case 3
        ratio=[]; colnum=4;
        for kk=unique(gIX)'
            subplot(ceil(max(unique(gIX))/colnum),colnum,kk,'ColorOrder',ColorOrder);
            if isplotback
                plot(mean(A(:,:,:),3),'linewidth',1,'color',[0.5 0.5 0.5]);hold on
            end
            %CS bar
            color=[0.5 0.5 0.5];
            patch1=patch(([frame.cs_start...
                frame.cs_end...
                frame.cs_end'...
                frame.cs_start]),...
                [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
                color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
            %US line
            l1=line(([frame.us_start frame.us_start]),[min(yli(:)) max(yli(:))],'color',[1 0 0],'linestyle','--','linewidth',2);hold on
            ind=cIX(find(gIX==kk)); ratio(kk,1)=length(ind);
            m=mean(A(:,:,ind),3);
            p1=plot(m,'linewidth',1.5);hold on
            ylim(yli);
            title([num2str(kk,'%02d')]);
            if kk==1
                legend([patch1,l1,p1'],{['CS';'US';num2str([1:size(m,2)]','%02d')]},'location','northwest'); 
            end
            trace_for_ttest{kk}=A(:,:,ind);
        end 
    case 4
                colnum=1;ratio=[];
        for kk=unique(gIX)'
            subplot(ceil(max(unique(gIX))/colnum),colnum,kk);
            if linevisible
                l3=plot(repmat([(frame.us_start):frame.per_cycle:size(A,1)],2,1)',[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',1.2);hold on
            end
            %CS bar
            color=[0.5 0.5 0.5];
            patch1=patch([[frame.cs_start:frame.per_cycle:size(A,1)]'...
                [frame.cs_end:frame.per_cycle:size(A,1)]'...
                [frame.cs_end:frame.per_cycle:size(A,1)]'...
                [frame.cs_start:frame.per_cycle:size(A,1)]']',...
                repmat([min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],size(A,1)/frame.per_cycle,1)',...
                color,'edgecolor',color,'facecolor',color,'edgealpha',0.4,'facealpha',0.45);hold on
            if isplotback
                plot(mean(reshape(A(:,:,:),size(A,1),[]),2),'LineStyle','--','linewidth',1.5,'color',[0.5 0.5 0.5]);hold on
            end
            ind=cIX(find(gIX==kk)); ratio(kk,1)=length(ind);
            m=squeeze(mean(A(:,:,ind),3));
            plot(m,'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            if length(ind)>1
            shadedErrorBar([1:length(m)],reshape(A(:,:,ind),size(A,1),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),size(A,1),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            end
            ylim(yli);title([num2str(kk,'%02d')]);
            set(gca,'xcolor','w','ycolor','w');
        end
end
set(groot,'defaultAxesColorOrder','remove');