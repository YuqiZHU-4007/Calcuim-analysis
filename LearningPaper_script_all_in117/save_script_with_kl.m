load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig_network\region_lab.mat');
 load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig_network\Plot_weight_selekon_sig_data\L_weight_selekton_net.mat')
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig_network\Plot_weight_selekon_sig_data\U_weight_selekton_net.mat')
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig_network\Plot_weight_selekon_sig_data\N_weight_selekton_net.mat')
group_label={'Learner (n=8)','Non-Learner (n=9)','Control (n=8)'};
sessionx = {'Pre Cond.','Early Cond,','Last Cond.','Post Cond.'};
savepath_eps='X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\';

A=cat(5,l_weight_selekton_net,N_weight_selekton_net,u_weight_selekton_net);
A=squeeze(A(:,1,:,:,:));s=size(A);A_r=reshape(A,1,[]);
A_n=zscore(A_r);A_n=reshape(A_n,s);
clr_cmapp =addcolorplus(302);
figure('position',[37,42,1775,953]),
ii=1;
for groupi=1:3
    for kk=1:4
        B=squeeze(A_n(kk,:,:,groupi));
        a=[-2 4];
        subplot(3,4,ii);imagesc(B,a);colormap(clr_cmapp);colorbar;
        set(gca,'XTick',[1:23],'XTickLabel',region_lab,'XTickLabelRotation',0);
        set(gca,'YTick',[1:23],'YTickLabel',region_lab,'YTickLabelRotation',0);
        set(gca, 'FontName', 'Arial', 'FontSize', 10);
        set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
        set(gcf,'Color','w')
        title(sessionx{kk});
        ii=ii+1;
    end
end


A=cat(5,l_weight_selekton_net,N_weight_selekton_net,u_weight_selekton_net);
clr_cmapp =addcolorplus(302);
ii=1;
for groupi=1:3
    B=squeeze(A(:,1,:,:,groupi));s=size(B);B_r=reshape(B,1,[]);
    B_n=zscore(B_r);B_n=reshape(B_n,s);
    for kk=1:4
        B=squeeze(B_n(kk,:,:));
        a=[-1 4];
        h=figure('position',[-1,42,633,555]);%subplot(3,4,ii);
        imagesc(B,a);colormap(clr_cmapp);colorbar;
        set(gca,'XTick',[1:23],'XTickLabel',region_lab,'XTickLabelRotation',90);
        set(gca,'YTick',[1:23],'YTickLabel',region_lab,'YTickLabelRotation',0);
        set(gca, 'FontName', 'Arial', 'FontSize', 12);
        set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
        set(gcf,'Color','w')
        title(sessionx{kk});
        name=[group_label{groupi},'_',sessionx{kk}];
        savefig(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.fig']);
        saveas(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.emf']);
        close(h)
    end
end



%%
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig\Plot_weight_selekon_sig_data\region_lab.mat');
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig\Plot_weight_selekon_sig_data\L_post_pre_sig.mat')
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig\Plot_weight_selekon_sig_data\U_post_pre_sig.mat')
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\Plot_fig\Plot_weight_selekon_sig_data\N_post_pre_sig.mat')
A=cat(3,l_post_pre_sig,n_post_pre_sig,u_post_pre_sig);
clr_cmapp =addcolorplus(301);
for groupi=1:3
    B=A(:,:,groupi)';
    h=bubbleheatmap(B,B,region_lab,region_lab,clr_cmapp,[0 1],[0 1.6]);colorbar('off');
    set(gca,'XTickLabelRotation',90);
    set(gca,'YTickLabelRotation',0);
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
    set(gcf,'Color','w')
    name=[group_label{groupi},'_','significant_map'];
    savefig(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.fig']);
    saveas(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.emf']);
    close(h)
end

aa=squeeze(sum(A,1))+squeeze(sum(A,2));

clr_cmapp =addcolorplus(301);
for groupi=1:3
    B=A(:,:,groupi);B=squeeze(sum(B,1))';
    h=bubbleheatmap(B,B,region_lab,region_lab,clr_cmapp2,[0 7],[0 1.6]);colorbar('off');
    set(gca,'XTickLabelRotation',90);
    set(gca,'YTickLabelRotation',0);
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
    set(gcf,'Color','w')
    name=[group_label{groupi},'_','significant_map'];
    savefig(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.fig']);
    saveas(h,[checkpath(fullfile(savepath_eps,'Network_with_KL')),'\',name,'.emf']);
    %close(h)
end

%%
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\网络图\Fish_post_pre_sig.mat')
h=figure('position',[-1,42,350,1800]);
for ii=1:4
A=squeeze(nanmean(n_fish_net(ii,:,:,:),4));B=squeeze(nanmean(n_fish_net(4,:,:,:),4));
clr_cmapp =addcolorplus(302);
a=[-1 4];
subplot(4,1,ii);
imagesc(A,[0 0.1]);colormap(clr_cmapp);colorbar;
set(gca,'XTick',[1:23],'XTickLabel',region_lab,'XTickLabelRotation',90);
set(gca,'YTick',[1:23],'YTickLabel',region_lab,'YTickLabelRotation',0);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
set(gcf,'Color','w')
end

p=[];h=[];
for ii=1:23
    for jj=1:23
        A=(squeeze(l_fish_net(1,ii,jj,:)));B=(squeeze(l_fish_net(4,ii,jj,:)));
        [p(ii,jj),h(ii,jj)]= ranksum(abs(A),abs(B));
    end
end
clr_cmapp =addcolorplus(301);
hh=bubbleheatmap(h,h,region_lab,region_lab,clr_cmapp,[0 1],[0 1.6]);colorbar('off');
set(gca,'XTickLabelRotation',90);
set(gca,'YTickLabelRotation',0);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
set(gca,'XColor','k','YColor','k','linewidth',2,'color','w');
set(gcf,'Color','w')
aa=squeeze(sum(h,1))+squeeze(sum(h,2))'; hh=bubbleheatmap(aa,aa,region_lab,region_lab,clr_cmapp,[0 8],[0 1.6]);colorbar('off');


