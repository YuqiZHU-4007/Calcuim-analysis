% ��������ϵ��
function [sil_values,silhouette_avg]= calculate_silhouette(data, eps, min_samples)
    [~, labels] = dbscan(data, eps,  min_samples,'Distance','precomputed');
    if any(labels == -1)
        % ����������㣬ֻ���Ǿ����
        core_samples =find( labels >= 0);
        labels = labels(core_samples);
        data = data(core_samples, :);
    end
    % ʹ��silhouette������������ϵ��
     sil_values = silhouette(data, labels);
     silhouette_avg = mean(sil_values(sil_values > -1));
end
