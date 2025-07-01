A = reshape(A(:,1:3075)',75,[],706);
B = smoothdata(A,1,'gaussian',5);
explained_var = sum(latent(1:2))/sum(latent);
[coeff,score,latent] = pca(squeeze(B(:,10,:))); 
neuro_subspace = coeff(:,1:3);

projected_data = zeros(41,75,3);
for i = 1:41
    
    projected_data(i,:,:) = squeeze(B(:,i,:))*neuro_subspace;
    plot3(squeeze(projected_data(i,:,1)),squeeze(projected_data(i,:,2)),squeeze(projected_data(i,:,3)));
    hold on
    
end