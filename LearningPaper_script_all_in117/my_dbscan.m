function [best_eps, best_min_samples, best_score,eps]=my_dbscan(data,dists,K,kD)

% 寻找合适的eps
minpts= K+1;
k = K^2-1; % 使用第2%的点
dists_sorted = sort(dists, k);
eps = mean(dists_sorted(:, round(0.02 * size(dists_sorted, 1))));
k_dist=sort(kD(end,:),'descend');
figure;plot(k_dist);hold on;
scatter(max(find(k_dist>=eps)),eps,'r','filled');
xlabel('点索引');
ylabel('第 k 近邻的距离');
title('第 k 近邻的距离图');
[idx_dbscan, corepts_dbscan] = dbscan(activities_dfdf_align_all_dist2,eps,minpts,'Distance','precomputed');
length(unique(idx_dbscan))
sum(idx_dbscan == -1)


% 寻找最优参数
eps_values = [eps];
min_samples_values = 2:20;
best_eps = 0;
best_min_samples = 0;
best_score = -1;
for i = 1:length(eps_values)
    eps = eps_values(i);
    for j = 1:length(min_samples_values)
        min_samples = min_samples_values(j);
        score = calculate_silhouette(data, eps, min_samples);
        if score > best_score
            best_score = score;
            best_eps = eps;
            best_min_samples = min_samples;
        end
    end
end
% 输出最优参数
fprintf('最优参数: eps=%.2f, min_samples=%d, 轮廓系数=%.2f\n', best_eps, best_min_samples, best_score);
