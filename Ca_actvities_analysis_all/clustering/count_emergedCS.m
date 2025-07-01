%% count of emerged CS responses
a=stimcorr.CSCS;normalize(stimcorr.CSCS,1,'zscore');
auto_thres=0.7;mean(a(:,1))+4*std(a(:,1));ind_emreges_CS={}; c=GetColormap('hsv_new',size(a,2));
figure,h=[];
for ii=1:size(a,2)
    h(ii,1)=histogram( a(:,ii),'BinEdges',-1:0.01:1,'FaceColor',c(ii,:));hold on
    %legend([h],{'Hab.','Tst.'});
    ind_emreges_CS{ii}=find(a(:,ii)>=auto_thres);
end
legend([h],{'Hab.','Acq.1','Acq.2','Acq.3','Acq.4','Acq.5','Tst.'});
line([auto_thres auto_thres],[0 1500],'linewidth',1.5,'color','r','linestyle','--');hold on
xlim([-1 1]);

clrmap_iii=[0 0 1;1 0 0];fraction_in_region_in_clust_per_block=nan(21,2,size(a,2),max(unique(index_all)'));
for ii=1:size(a,2)
    cIX_iii=1:size(a,1);gIX_iii=ones(size(a,1),1);gIX_iii(ind_emreges_CS{ii})=2;zz=1;
    for kk=unique(index_all)'
        ind_in_this_fish=(find(index_all==kk));
        ind_in_this_clus=find(gIX_iii== 2);%1:length(gIX_iii);%
        cIX_in_this_clus=cIX_iii(ind_in_this_clus);%
        [IX,i_ind_in_this_fish,i_cIX]=intersect(ind_in_this_fish' ,cIX_in_this_clus');
        %unique(index_all(IX))
        %     unique(index_all(cIX_iii(ind_in_this_clus(i_cIX))))
        %     unique(gIX_iii(ind_in_this_clus(i_cIX)))
        %e=load(envbatch{kk});
        load([savebatch{kk} '\brain arearegion_mask.mat']);
        [~,b,~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),cIX_iii(ind_in_this_clus(i_cIX)),...
            [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
            [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap_iii,false);
        
        title(['Fish ' num2str(kk,'%02d')]);
        if ~isempty(b)
            fraction_in_region_in_clust_per_block(:,:,ii,zz)=b;
        end
        zz=zz+1;
    end
end

%plot farction of emergeCS per session
figure,
Lable={'l pallium rostral'; 'r pallium rostral';'l pallium middle';'r pallium middle';'l pallium lateral';'r pallium lateral';...
    'l habenula';'r habenula';...
    'l tectum rostral';'r tectum rostral';'l tectum middle'; 'r tectum middle';'l tectum caudal';'r tectum caudal';...
    'l tegmentum';'r tegmentum';...
    'l cerebellum';'r cerebellum';...
    'l hindbrain';'r hindbrain';...
    'longitude'
    };c = categorical(Lable)';
for tt=1:size(a,2)
    subplot(size(a,2),1,tt),
    l_fraction_in_region_in_clust=mean(fraction_in_region_in_clust_per_block(:,2,tt,1:4),4,'omitnan');
    nl_fraction_in_region_in_clust=mean(fraction_in_region_in_clust_per_block(:,2,tt,5:end),4,'omitnan');
    d_fraction_in_region_in_clust=[l_fraction_in_region_in_clust,nl_fraction_in_region_in_clust];
    
    if tt<size(a,2)
        h=bar(d_fraction_in_region_in_clust,'FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);set(gca,'xticklabel',[]);
    else
        h=bar(c,d_fraction_in_region_in_clust,'FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);
        set(gca,'xticklabels',Lable,'XTickLabelRotation',45);
        legend('Learner','Non-learner');
    end
    h(1).CData = [1 0 0];h(2).CData = [0 0 1];
    title([num2str(tt,'%02d')]);ylim([0 0.2])
end
