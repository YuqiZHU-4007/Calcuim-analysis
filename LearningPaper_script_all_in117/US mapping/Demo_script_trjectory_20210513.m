%% trojectory
%trojectory of newly emerged
clc
clear all;

%correlation map
n=10; figure;
for ii=size(activities_spon,3):-1:1
    a=squeeze(dfdf(CS_new_emerged,:,ii));
    Cor=corr(a');
    if ii==size(activities_spon,3)
        %figure,imagesc(Cor,[0 max(Cor(:))]);colormap('hot');colorbar;
        [idx_functional,~,~,~] = kmeans(Cor',n,'Distance','correlation','Replicates',10);
        [idxAs,idxAa_ind]=sort(idx_functional);
    end
    %figure;hist(idx_functional,[1:n]);xlim([1 n]);
    subplot(2,5,ii);imagesc(Cor(idxAa_ind,idxAa_ind),[-0.2 0.2]);colormap('jet');colorbar;%axis equal
    title(num2str(ii))
end
indd=CS_new_emerged(idxAa_ind);
figure,scatter3(supervoxel(indd,1),supervoxel(indd,2),supervoxel(indd,3),2,[0.5 0.5 0.5],'filled');axis equal;
%PCA
A1=zeros(size(dfdf,2),1,length(CS_new_emerged));A2=zeros(size(dfdf,2),1,length(CS_new_emerged));
A1(:,1,:)=dfdf(CS_new_emerged,:,1)';%A1 = reshape(A1',frame.per_cycle,[],length(CS_new_emerged));
%A1=mean(A1,2);
A2(:,1,:)=dfdf(CS_new_emerged,:,end)';%A2 = reshape(A2',frame.per_cycle,[],length(CS_new_emerged));
%A2=mean(A2,2);
[projected_data,neuro_subspace,explained_var]=mypca(cat(2,A2,A1),20,0,frame);

%t-sne