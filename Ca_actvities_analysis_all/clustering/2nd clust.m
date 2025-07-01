num_clust=[1:41];ind_in_num_clust=[];zz=1;a=[];
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
yli=[-0.02 0.04];kk=1;pt=[];ht=[];ind_in_num_clust=[];
 for jj=num_clust
    ind_in_num_clust=cat(1,ind_in_num_clust,find(gIX_iii==jj));
    %correct 2
%     if jj==4 | jj==7
%         [ind_c,~,i_ind]=intersect(find(h_acq(:,1)==1),cIX_iii(ind));ind=ind(i_ind);
% % %     elseif jj==10
% % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==1),cIX_iii(ind));ind=ind(i_ind);
% % %     elseif jj==20
% % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==1),cIX_iii(ind));ind=ind(i_ind);
% % %     elseif jj==21
% % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,3)==1),cIX_iii(ind));ind=ind(i_ind);
%     elseif jj==23
%         [ind_c,~,i_ind]=intersect(find(h_acq(:,1)==0),cIX_iii(ind));ind=ind(i_ind);
%         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==0),cIX_iii(ind));ind=ind(i_ind);
%         [ind_c,~,i_ind]=intersect(find(h_acq(:,5)==1),cIX_iii(ind));ind=ind(i_ind);
%      end
    a(zz)=length(ind_in_num_clust);  zz=zz+1;
 end
rank_ind=find(h_base_adj==1 & h_tst_adj==0);rank_ind=find(h_tst_adj==1);rank_ind=find(h_base_adj==0 & h_tst_adj==1);%i_ind=i_ind(500);
[ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(ind_in_num_clust));ind_in_num_clust=ind_in_num_clust(i_ind);
%cur motor related
%figure,hist(motorcorr.sep_CS(cIX_iii(ind),7))
ind_in_num_clust(find(motorcorr.sep_CS(cIX_iii(ind_in_num_clust),7)>=0.3 ))=[];%& motorcorr.sep_spon(cIX_iii(ind),7)>=0.4
figure,plot(mean(squeeze(act_all_hab_mean(:,:,cIX_iii(ind_in_num_clust))),2));hold on;plot(mean(squeeze(act_all_tst_mean(:,:,cIX_iii(ind_in_num_clust))),2))
c=7;
figure,
for ii=1:trial.acq_block_num+2
    subplot(1,trial.acq_block_num+2,ii),
    if ii==1
        a=reshape(act_all_hab_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),1,[]);b=squeeze(a(:,1,cIX_iii(ind_in_num_clus)))';
    elseif ii==7
        a=reshape(act_all_tst_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),1,[]);b=squeeze(a(:,1,cIX_iii(ind_in_num_clus)))';
    else
        a=reshape(act_all_acq_block_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*trial.acq_block_num,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),trial.acq_block_num,[]);b=squeeze(a(:,ii-1,cIX_iii(ind_in_num_clus)))';
        %a=act_all_acq_block_mean(frame.cs_start-c:frame.us_start-1,:,:);%
    end
    b=smoothdata(b,1,'gaussian' ,5);
    imagesc(b);hold on;
    line([c c]+0.5,[0 length(ind_in_num_clus)],'color','r','linestyle','--','linewidth',2);
    set(gca,'clim',[0.5 1.5]);%set(gca,'visible','off')
end

for ii=1
    %% autocluster
    addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
    frame_ind=frame.cs_start:frame.us_start-2;%%%%%%%%%%%%%%%%%%%!!!!!!!
    M_0=reshape(act_all_acq_block_mean(frame_ind,:,cIX_iii(ind_in_num_clust)),length(frame_ind)*size(act_all_acq_block_mean,2),[])';
    M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
    %M_0=regressor_profile_n_p(cIX_iii(ind_in_num_clust),:);
    %%%clustering
    numK=20;
    %numK=max(idx_kmeans_seed);
    % ind=index_all;%find(index_all==2);
    % ind=ind([1:20:length(ind)]);
    % [gIXx,~,numK,~]=find_best_k_in_range(M_0(ind,:),3:15);
    % length(unique(gIXx))
    % cIXx=ind;cIX_reg = (1:size(M_0,1))';%ind;
    if isempty(numK)
        numK=20;
    end
    masterthres=0.6;
    para=struct;
    para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
    %cIX=(1:size(M_0,1))';gIX=(1:size(M_0,1))'; cIX_reg = (1:size(M_0,1))';
    cIXx_2nd=1:size(M_0,1);gIXx_2nd=gIX_iii(ind_in_num_clust);cIX_reg_2nd = (1:size(M_0,1))';
    for tt=1:5
        %[cIX,gIX] = AutoClustering(randi(size(M_0,1),[1,50000]),[1:50000],M_0,cIX_reg,1,para,1,0.7);
        [cIX_2nd,gIX_2nd] = AutoClustering(cIXx_2nd,gIXx_2nd,M_0,cIX_reg_2nd,0,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
        length(unique(gIX_2nd))
        ind=[1:length(cIX_2nd)];%randi(length(cIX),[floor(length(cIX)/5),1]);
        cIXx_2nd=cIX_2nd(ind);
        gIXx_2nd=gIX_2nd(ind);
    end
    length(unique(gIX_2nd))
    type_2ndclus.ind=ind_in_num_clust;
    type_2ndclus.cIXx_2nd=cIXx_2nd;
    type_2ndclus.gIXx_2nd=gIXx_2nd;
end
gIX_2nd_to_all=type_2ndclus.gIXx_2nd;
cIX_2nd_to_all=cIX_2nd_to_all(ind_in_num_clust(type_2ndclus.cIXx_2nd));
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX));
ind_control=setdiff(1:length(index_all),cIX_2nd_to_all)';%length(ind)+length(cIX)== length(index_all) 
vi='on';'off';
for iii=1
    switch iii
        case 1
            cIX_2nd_to_all=cIX_2nd_to_all;gIX_2nd_to_all=gIX;clrmap_iii=clrmap;
    end
    outputpath=checkpath([savepath_all '\figures_regressor_profile' '_' num2str(iii,'%02d')]);
    %count
    figure,h=histogram(gIX_2nd_to_all,'binedge',[1:max(unique(gIX_2nd_to_all))+1]);h=bar(h.Values/length(gIX_all));h.FaceColor = 'flat';
    for ii=unique(gIX_2nd_to_all)
        h.CData(ii,:) = clrmap_iii(ii,:);
    end
    xlabel('# Cluster');ylabel('Fraction');
    %1
    [h,gIX_2, numU] = hierplot_zyq_20190530(cIX_2nd_to_all,gIX_2nd_to_all,regressor_profile_n_p(cIX_2nd_to_all,:));saveas(h,[outputpath '\f1'],'fig');
    figure,imagesc(regressor_profile_n_p(cIX_2nd_to_all,:) ,[-2 2]);colorbar;colormap('parula');
    %2
    a=gIX_2nd_to_all;a(find(gIX_2nd_to_all==13))=[];b=cIX_2nd_to_all;b(find(gIX_2nd_to_all==13))=[];
     B=act_all_acq(frame.cs_start:frame.us_start-1,:,:,:);B=reshape(B,length([frame.cs_start:frame.us_start-1])*trial.acq(1),[])';
     [h]=pushbutton_popupplot_Callback(B,b,a,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[outputpath '\f2'],'fig');
     B=act_all_acq_block_mean(frame.cs_start:frame.us_start-2,:,:,:);B=reshape(B,length([frame.cs_start:frame.us_start-2])*trial.acq_block_num,[])';
     [h]=pushbutton_popupplot_Callback(B,b,a,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);saveas(h,[outputpath '\f2'],'fig');
    %3
%     h=figure;B=cat(1,area_all.CS_hab_tst(1,:),area_all.CS_acq_block,area_all.CS_hab_tst(2,:));%B=[area.CS_hab_tst(1,:); area.CS_acq_block ;area.CS_hab_tst(2,:)];
%     a=[-0.01 1];plot_test_1([1:size(B,1)],B,cIX_iii,gIX_iii,clrmap_iii,a,frame,fs.ca,2);saveas(h,[outputpath '\f3'],'fig');
    %4
    %[h,ratio,~]=plot_test_2(act_all_acq_block_mean,cIX_iii,gIX_iii,frame,clrmap_iii,1,[-0.02 0.04],true,true);saveas(h,[outputpath '\f41'],'fig');
    [h,ratio,~]=plot_test_2(act_all_acq_block_mean,cIX_2nd_to_all,gIX_2nd_to_all,frame,clrmap_iii,2,[-0.02 0.04],true,true);saveas(h,[outputpath '\f42'],'fig');
    [h,ratio,trace_for_ttest_acq]=plot_test_2(act_all_acq_block_mean,cIX_2nd_to_all,gIX_2nd_to_all,frame,clrmap_iii,3,[-0.005 0.02],true,true);saveas(h,[outputpath '\f42'],'fig');
    %5
    [h,~,trace_for_ttest_hab]=plot_test_2(act_all_hab_mean,cIX_2nd_to_all,gIX_2nd_to_all,frame,clrmap_iii,2,[-0.01 0.02],false,true);saveas(h,[outputpath '\f5'],'fig');
    %6
    [h,~,trace_for_ttest_test]=plot_test_2(act_all_tst_mean,cIX_2nd_to_all,gIX_2nd_to_all,frame,clrmap_iii,2,[-0.01 0.02],false,true);saveas(h,[outputpath '\f6'],'fig');
    %
    plot_test_4(act_all_hab_mean,act_all_tst_mean,cIX_2nd_to_all,gIX_2nd_to_all,frame,clrmap_iii,2,[-0.01 0.02],false,true);
    %7
    [h1,h2,ratio_fish]=plot_test_3(act_all',cIX_2nd_to_all,gIX_2nd_to_all,clrmap_iii,index_all,envbatch,env_all,outputpath,fs,stimCS,stimUS);
end
