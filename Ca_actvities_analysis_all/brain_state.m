%% brain state
%% PCA
bin=1:15;
A=graph.node_mean_act;
A=smoothdata(A',2,'movmean',20);
A=A(1:3075,:);
A=reshape(A,frame.per_cycle,[],706);
A=squeeze(mean(A,2));
%B=mean(A,1);
figure,imagesc(A',[0 0.05]);colorbar;
%A=normalize(A,1,'zscore');
[cof,pca1,latent,tsquare]  = pca(A);
figure,plot(cumsum(latent)./sum(latent))

A=graph.node_mean_act2;
A=smoothdata(A',2);
figure,imagesc(A',[0 0.05]);colorbar;
A=normalize(A',2,'zscore');
[cof,pca2,latent,tsquare]  = pca(A');
figure,plot(cumsum(latent)./sum(latent))

figure,plot3(pca1(:,1),pca1(:,2),pca1(:,3),'r');hold on;
plot3(pca2(:,1),pca2(:,2),pca2(:,3),'b');hold on;
xlabel('pc1'); ylabel('pc2'); zlabel('pc3');
%%
% spon1 bef VS aft
A1=graph.node_mean_act(:,1:3000);A1 = reshape(A1',frame.per_cycle,[],size(graph.nodes,1));A1=mean(A1,2);
A2=graph.node_mean_act2(:,1:3000);A2 = reshape(A2',frame.per_cycle,[],size(graph.nodes,1));A2=mean(A2,2);
%A3=graph.node_mean_act_training;A3 = reshape(A3',frame.per_cycle,[],size(graph.nodes,1));
A=cat(2,A2,A1);
[projected_data,neuro_subspace,explained_var]=mypca(cat(2,A2,A1),0);
c=[1 0 0; 0 0 1];
% figure;
% for i = 1:size(A,2)
%     if i<40
%     plot3(squeeze(projected_data(i,:,1)),squeeze(projected_data(i,:,2)),squeeze(projected_data(i,:,3)),'color',c(1,:));hold on;
%     else
%         plot3(squeeze(projected_data(i,:,1)),squeeze(projected_data(i,:,2)),squeeze(projected_data(i,:,3)),'color',c(2,:));hold on;
%     end
%     %plot3(squeeze(score(i,:,1)),squeeze(score(i,:,2)),squeeze(score(i,:,3)));hold on
% end
% xlabel('pc1'); ylabel('pc2'); zlabel('pc3');
%
%% tsne
M=normalize(graph.node_mean_act_training,2,'zscore');
corr=re_startpoint(find(re_startpoint(:,2)>=frameb.cs_start & re_startpoint(:,2)<= frameb.cs_end ),1);
%re_startpoint(:,1)>= trial.test(2) & re_startpoint(:,1) <= trial.test(1)
corr=intersect(corr,[trial.test(2):trial.test(3)])';
incorr=setdiff([trial.test(2):trial.test(3)],corr);
a=40;
for  perplexity=a
%     figure;
    numDims = 2; mappedA_corr=[];
    for ii=corr
        ind=(ii-1)*frame.per_cycle+1:ii*frame.per_cycle;
        Mm=M(:,ind);
        %mappedA_corr(:,:,ii) = tsne(Mm', [], numDims,[],perplexity ); % no_dims = 2
        [~,mappedA_corr(:,:,ii),~,~]  = pca(Mm');
%         plot(mappedA_corr(:,1,ii),mappedA_corr(:,2,ii),'linewidth',1);hold on;
    end
%     title(num2str(perplexity));
end
for  perplexity=a
%     figure;
    numDims = 2; mappedA_incorr=[];
    for ii=incorr
        ind=(ii-1)*frame.per_cycle+1:ii*frame.per_cycle;
        Mm=M(:,ind);
%         mappedA_incorr(:,:,ii) = tsne(Mm', [], numDims,[],perplexity ); % no_dims = 2
     [~,mappedA_incorr(:,:,ii),~,~]  = pca(Mm');
%         plot(mappedA_incorr(:,1,ii),mappedA_incorr(:,2,ii),'linewidth',1);hold on;
    end
%     title(num2str(perplexity));
end
 mappedA_corr(:,:,find( mappedA_corr(1,1,:)==0))=[];mappedA_corr_mean=mean(mappedA_corr,3);
 mappedA_incorr(:,:,find( mappedA_incorr(1,1,:)==0))=[];mappedA_incorr_mean=mean(mappedA_incorr,3);
figure;
for ii=1:size(mappedA_corr,3)
    plot3(mappedA_corr(:,1,ii),mappedA_corr(:,2,ii),mappedA_corr(:,3,ii),'color',[0.5 0.5 0.5],'linewidth',1);hold on;
end
plot3(mappedA_corr_mean(:,1),mappedA_corr_mean(:,2),mappedA_corr(:,3,ii),'color',[1 0.1 0.1],'linewidth',2);hold on;
for ii=1:size(mappedA_incorr,3)
    plot3(mappedA_incorr(:,1,ii),mappedA_incorr(:,2,ii),mappedA_incorr(:,3,ii),'color',[0.5 0.5 0.5],'linewidth',1);hold on;
end
plot3(mappedA_incorr_mean(:,1),mappedA_incorr_mean(:,2),mappedA_incorr(:,3,ii),'color',[0.1 0.1 1],'linewidth',2);hold off;