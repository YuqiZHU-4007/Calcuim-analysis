function [h,AUC,AUC_go,AUC_nogo,p_AUC,p_AUC_go,p_AUC_nogo,CR_ratio]=plot_AUC(act_all,CR,trial,frame)
CR_ratio=[];AUC=[];AUC_go=[];AUC_nogo=[];
frame_ind_cs=frame.cs_start:frame.cs_end-1;
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            trial_ind=trial.hab(2):trial.hab(3);
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            trial_ind=(trial.hab(3)+(ii-2)*trial.acq_block_trial)+1 :(trial.hab(3)+(ii-1)*trial.acq_block_trial);
        case trial.acq_block_num+2
            trial_ind=trial.test(2):trial.test(3);
    end
    
    go_trial=unique(CR(find(CR(:,1)<=trial_ind(end) & CR(:,1)>=trial_ind(1)),1));
    nogo_trial=setdiff(trial_ind,go_trial);
    CR_ratio(ii,1)=length(go_trial)./length(trial_ind);
    AUC(ii,:)=squeeze(mean(sum(act_all(frame_ind_cs,trial_ind,:),1,'omitnan'),2,'omitnan'));
    AUC_go(ii,:)=squeeze(mean(sum(act_all(frame_ind_cs,go_trial,:),1,'omitnan'),2,'omitnan'));
    AUC_nogo(ii,:)=squeeze(mean(sum(act_all(frame_ind_cs,nogo_trial,:),1,'omitnan'),2,'omitnan'));
end
p_AUC=[];p_AUC_go=[];p_AUC_nogo=[];
for ii=1:trial.acq_block_num+2
    a=AUC(1,:);b=AUC(ii,:);
    if sum(isnan(b))~=length(b)
        p_AUC(ii,:)=ranksum(a,b);
    else
        p_AUC(ii,:)=nan;
    end
    
    a=AUC_go(1,:);b=AUC_go(ii,:);
    if sum(isnan(a))==length(a) a=AUC(1,:);end
    if sum(isnan(b))~=length(b)
        p_AUC_go(ii,:)=ranksum(a,b);
    else
        p_AUC_go(ii,:)=nan;
    end
    
    a=AUC_nogo(1,:);b=AUC_nogo(ii,:);
    if sum(isnan(a))==length(a) a=AUC(1,:);end
    if sum(isnan(b))~=length(b)
        p_AUC_nogo(ii,:)=ranksum(a,b);
    else
        p_AUC_nogo(ii,:)=nan;
    end
end

h=figure('position',[124,558,1559,420]);t=tiledlayout(1,3);%t.TileSpacing = 'compact';%t.Padding = 'compact';
for ii=1:3
    ax1=nexttile([1,1]);
    switch ii
        case 1
            AAA=AUC';CCC=[1 0 0];p=p_AUC;seg='All trials';
        case 2
            AAA=AUC_go';CCC=[1.00,0.41,0.16];p=p_AUC_go;seg='Go trials';
        case 3
            AAA=AUC_nogo';CCC=[1.00,0.07,0.65];p=p_AUC_nogo;seg='No-Go trials';
    end
    yyaxis left;shadedErrorBar([1:trial.acq_block_num+2],AAA,{@(x) mean(x,'omitnan'),@(x) std(x)*(1/sqrt(size(AAA,1)))},'lineprops',{CCC},'transparent',1,'patchSaturation',0.1); hold on
    ylabel('AUC','fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
    yyaxis right;
    plot([1:trial.acq_block_num+2], CR_ratio,'k-.','linewidth',1.5,'markersize',8);hold on
    ylabel('CR ratio','fontsize',16,'FontWeight','bold');ylim([0 1.2]);set(gca,'ycolor',[0 0 0],'linewidth',2);
    a=find(p<0.05); yyaxis right;scatter(a,0.05*ones(length(a),1),16,'k','*');
    xlim([1 trial.acq_block_num+2]);set(gca,'xtick',[1:trial.acq_block_num+2],'xticklabel',{'Pre','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post'},'XTickLabelRotation',45,'fontsize',12);
    title(seg,'fontsize',16,'FontWeight','bold');
end
