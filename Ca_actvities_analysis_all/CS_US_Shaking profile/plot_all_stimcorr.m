clc;clear all;
load('K:\3.Poster_fig\part3\non-learner\stAvrcorr_all_raw.mat')
savepath_eps='K:\3.Poster_fig\part3\non-learner\';
savetype='epsc';
stAvrcorr_all=[];
n=4;
for ii=3%n:size(stAvrcorr_all_raw,1)
    stAvrcorr_all=[stAvrcorr_all,stAvrcorr_all_raw{ii,1}];
end
h5=figure('name',['Num.','_hist_of_StimAveCorr_autothr','_ALL']);
c=GetColormap('hsv_new',5);%[0 0 1;1 0 0,];
kk=1;h=[];
for ii=[2:6]
    [~,IX] = sort(stAvrcorr_all(ii,:)','descend');
    h(kk,1)=histogram( stAvrcorr_all(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:),'Normalization','probability');hold on
    %line([auto_thres_bef_conditioing auto_thres_bef_conditioing],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
    kk=kk+1;
end
legend([h],{'Hab.','Acq.1','Acq.2','Acq.3','Acq.4','Acq.5','Tst.'});
%legend([h],{'Hab.','Tst.'});
xlim([-1 1]);set(gca,'fontsize',20)
saveas(h5,[savepath_eps 'Num.','_hist_of_StimAveCorr_autothr','_ALL' '.eps'],savetype);
%%
h6=figure('name',['Num.','_histcount_of_StimAveCorr_autothr','_ALL']);
c=GetColormap('hsv_new',7);%[0 0 1;1 0 0,];
c=[0.93,0.38,0.28;0.91,0.91,0.07;0.06,1.00,1.00;0.06,1.00,1.00;1.00,0.07,0.65;1 0 0];
%c=[0 0 1;1 0 0];
kk=1;h=[];
for ii=[2:6]
    [~,IX] = sort(stAvrcorr_all(ii,:)','descend');
    hh=histogram( stAvrcorr_all(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:),'Normalization','cdf','visible','off');
    y=hh.Values;x=hh.BinEdges(2:end);
    h(kk,1)=plot(x,y,'color',c(kk,:),'linewidth',2);hold on;
    %line([auto_thres_bef_conditioing auto_thres_bef_conditioing],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
    kk=kk+1;
end
%legend([h],{'Hab.','Acq.1','Acq.2','Acq.3','Acq.4','Acq.5','Tst.'});
xlim([-1 1]);ylim([0 1]);
set(gca,'fontsize',20)
saveas(h6,[savepath_eps 'Num.','_histcount_of_StimAveCorr_autothr','_ALL' '.eps'],savetype);

%%
c=GetColormap('hsv_new',7);%[0 0 1;1 0 0,];
x=[];y=[];
for jj=n:size(stAvrcorr_all_raw,1)
    kk=1; X=stAvrcorr_all_raw{jj,1};
    for ii=[1:7]
        [~,IX] = sort(X(ii,:)','descend');
        hh=histogram( X(ii,:),'BinEdges',-1:0.01:1,'FaceColor',c(kk,:),'Normalization','cumcount','visible','off');
        y(ii,:,jj)=hh.Values;x(ii,:,jj)=hh.BinEdges(2:end);
        kk=kk+1;
    end
end
h7=figure('name',['Num.','_histcount_of_StimAveCorr_autothr','_ALL_errorbar']);
kk=1;h=[]; xx=mean(x,3);yy=mean(y,3);errBar=std(y,0,3)/size(y,3);
for ii=[1:7]
    %errBar=std(reshape(y(ii,:,:),[],size(stAvrcorr_all_raw,1)),0,3);
    shadedErrorBar(xx(ii,:),yy(ii,:),errBar(ii,:),'lineProps',{c(kk,:)});hold on;
    %line([auto_thres_bef_conditioing auto_thres_bef_conditioing],[0 1500],'linewidth',1.5,'color',c(kk,:),'linestyle','--');hold on
    kk=kk+1;
end
xlim([-1 1]);
saveas(h7,[savepath_eps 'Num.','_histcount_of_StimAveCorr_autothr','_ALL_errorbar' '.eps'],savetype);
x=[];
for jj=n:size(stAvrcorr_all_raw,1)
X=stAvrcorr_all_raw{jj,1};
x(:,jj)=mean(X,2);
end
h8=figure('name',['Num.','_boxplot_of_StimAveCorr_autothr','_ALL']);
boxplot(flip(stAvrcorr_all',2),'Labels',{'Tst.','Acq.5','Acq.4','Acq.3','Acq.2','Acq.1','Hab.'},'Symbol','wo','OutlierSize',2,'Orientation','horizontal');
xlim([-0.4 0.8]);
saveas(h8,[savepath_eps,'Num.','_boxplot_of_StimAveCorr_autothr','_ALL' '.eps'],savetype);

h9=figure('name',['Num.','_boxplot_of_StimAveCorr_autothr','_ALL_mean']);
boxplot(flip(x',2),'Labels',{'Tst.','Acq.5','Acq.4','Acq.3','Acq.2','Acq.1','Hab.'},'Symbol','wo','OutlierSize',2,'Orientation','horizontal');
xlim([-0.1 0.5]);
saveas(h9,[savepath_eps,'Num.','_boxplot_of_StimAveCorr_autothr','_ALL_mean' '.eps'],savetype);

