function [h,ratio,trace_for_ttest]=plot_test_4(A,B,cIX,gIX,frame,clrmap,type,yli,linevisible,isplotback)
if  ~exist('type','var')
    type=1;
end
if length(size(A))<3
    A=reshape(A,[size(A,1),1,size(A,2)]);
    B=reshape(B,[size(B,1),1,size(B,2)]);
end
trace_for_ttest={};
ColorOrder=[0,0,0;...
    0,1,0;...
    1,0,0;...
    0,0,1;...
    1,1,0;];
set(groot,'defaultAxesColorOrder',ColorOrder)
h=figure;
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
            plot(m,'LineStyle' ,'--','linewidth',1.2,'color',clrmap(kk,:,:));hold on
            plot(mean(reshape(B(:,:,ind),frame.per_cycle*size(B,2),[]),2),'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1); hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1); hold on
            %plot(reshape(A(:,:,ind),frame.per_cycle*size(A,2),[]),'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            %plot_rawtrace(m,[1:length(m)],[],fs,frame,trial,[]);hold on;
        end
        ylim(yli);
    case 2
        colnum=4;ratio=[];
        for kk=unique(gIX)'
            subplot(ceil(max(unique(gIX))/colnum),colnum,kk);
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
            if isplotback
                plot(mean(reshape(A(:,:,:),frame.per_cycle*size(A,2),[]),2),'linewidth',1,'color',[0.5 0.5 0.5]);hold on
            end
            ind=cIX(find(gIX==kk)); ratio(kk,1)=length(ind);
            m=mean(reshape(A(:,:,ind),frame.per_cycle*size(A,2),[]),2);
            plot(m,'LineStyle' ,'--','linewidth',1.2,'color',clrmap(kk,:,:));hold on
            plot(mean(reshape(B(:,:,ind),frame.per_cycle*size(B,2),[]),2),'linewidth',1.2,'color',clrmap(kk,:,:));hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            %shadedErrorBar([1:length(m)],reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])',{@mean,@(x) std(x)*(1/sqrt(length(ind)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            trace_for_ttest{kk}=reshape(A(:,:,ind),frame.per_cycle*size(A,2),[])';
            %yli=[0 0.1];
            ylim(yli);title([num2str(kk,'%02d')]);
        end
end
set(groot,'defaultAxesColorOrder','remove');