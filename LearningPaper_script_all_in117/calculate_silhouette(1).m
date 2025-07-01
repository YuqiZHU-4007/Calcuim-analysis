% 计算轮廓系数
function [sil_values,silhouette_avg]= calculate_silhouette(data, eps, min_samples)
    [~, labels] = dbscan(data, eps,  min_samples,'Distance','precomputed');
    if any(labels == -1)
        % 如果有噪声点，只考虑聚类点
        core_samples =find( labels >= 0);
        labels = labels(core_samples);
        data = data(core_samples, :);
    end
    % 使用silhouette函数计算轮廓系数
     sil_values = silhouette(data, labels);
     silhouette_avg = mean(sil_values(sil_values > -1));
end
