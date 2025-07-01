function [cluster_ind_cut,cluster_ind_cut_remat,ind_cut,ind_save,value_cut,value_save]=cut_small_cluster(cluster_ind,n_cluster)
t=tabulate(cluster_ind);
value_cut=t(find(t(:,2)<n_cluster),1);value_save=t(find(t(:,2)>=n_cluster),1);
[ind_cut,~] = find(bsxfun(@eq,cluster_ind,value_cut'));
ind_save=setdiff(1:length(cluster_ind),ind_cut);
cluster_ind_cut=cluster_ind;cluster_ind_cut(ind_cut)=[];
cluster_ind_cut_remat=cluster_ind_cut;kk=1;
for ii=value_save'
    cluster_ind_cut_remat(find(cluster_ind_cut==ii))=kk;
    kk=kk+1;
end

