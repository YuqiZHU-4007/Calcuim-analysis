function [stimcorr_k,stAvrcorr_k,IX_passX]=get_regressor_M_corr(M,regress,opt,frame,frame_ind_cs)
prct_const=opt.prct_const;
nCells_total = size(M,1);
topN = round(prct_const/100 * nCells_total);
reg_thres=opt.reg_thres;reg_thres2=opt.reg_thres2;%para
colorCS=[0.5 0.5 0.8];
ii=1;
[stimcorr_k(:,ii),~] = MotorSourceCorrelation(M,regress,[]);
[~,IX] = sort(stimcorr_k(:,ii),'descend');
thr = stimcorr_k((IX(topN)),ii);
IX_passX = (find(stimcorr_k(:,ii)>max(reg_thres,thr)))';
stim_output=mean(M((IX_passX),:),1);
[stim_tAvr,~,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind_cs));
stAvrcorr_k(:,ii) = corr(stim_tAvr',M');
% for pp=1
%     [~,IX] = sort(stAvrcorr_k(:,ii),'descend');
%     x0 = stAvrcorr_k((IX(topN)),ii);
%     thr=max(x0,reg_thres2);IX_passX = find( stAvrcorr_k(:,ii)>=thr);
%     figure('name',['regressor vs comt temp ' num2str(ii)]),
%     histogram(stimcorr_k(:,ii),'BinEdges',-1:0.01:1,'FaceColor','b');hold on
%     histogram(stAvrcorr_k(:,ii),'BinEdges',-1:0.01:1,'FaceColor','r');hold on
%     figure('name',num2str(ii)),
%     x=stim_output';a=[-0.05 0.4];
%     patch1=patch([[frame.cs_start:length(frame_ind_cs):size(x,1)]'...
%         [min(frame.cs_end,frame_ind_cs(end)):length(frame_ind_cs):size(x,1)]'...
%         [min(frame.cs_end,frame_ind_cs(end)):length(frame_ind_cs):size(x,1)]'...
%         [frame.cs_start:length(frame_ind_cs):size(x,1)]']',...
%         repmat([min(a) min(a) max(a) max(a)],length([frame.cs_start:length(frame_ind_cs):size(x,1)]),1)',...
%         colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
%     p1=plot(reshape(stim_output,1,[]),'k','linewidth',1);hold on;
%     p2=plot(stim_tAvr,'b','linewidth',1);hold on;
%     %     if isus
%     %         l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
%     %     end
%     p3=plot(mean(M(IX_passX,:),1),'r','linewidth',1.5);hold on;
%     legend([patch1 p1,p2,p3],{'CS','Regressor','Comp. Templete','>thr Avg. trace'},'location','northeastoutside');
%     ylim(a);xlim([1 length(x)]);
%     box off;
%     %line([thr thr],[0 1500],'linewidth',1.5,'color','r','linestyle','--');hold on
% end