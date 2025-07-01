function [Mahalanobis_d_raw, Mahalanobis_mean]=my_mahal(reference_data, target_data)

% reference_data=UR_ds_PV(:,ind);%t * n
% target_data=CR_ds_PV(:,ind);%t * n
Mahalanobis_d_raw={};Mahalanobis_mean=nan(1,4);


%kaoru_2022_1
try
    [Mahalanobis_d_raw{1},  Mahalanobis_mean(:,1)] = fxn_nmt_Mahalanobis(reference_data, target_data);
catch
    warning('wrong');
end
%kaoru_2022_2
try
    [Mahalanobis_d_raw{2},  Mahalanobis_mean(:,2)] =  fxn_nmt_Mahalanobis(mean(reference_data,1,'omitnan'), mean(target_data,1,'omitnan'));
catch
    warning('wrong');
end

% %schnitzer_2017_3
% Mahalanobis_d_raw{4} = mahal(target_data,reference_data);% mean of refer data
% Mahalanobis_mean{4} = mean(Mahalanobis_d_raw{4});
%schnitzer_2017_1
try
    Mahalanobis_d_raw{3} = mahal(mean(target_data,1,'omitnan')',mean(reference_data,1,'omitnan')');% mean of refer data
    Mahalanobis_mean(:,3) = mean(Mahalanobis_d_raw{3});
catch
    warning('wrong');
end
%kaoru_2022_2
try
    [neuro_subspace,~, target_data_PCA_score] = fxn_Marchenko_thrcover_PCA(target_data, 90);
    reference_data_PCA_score = reference_data*neuro_subspace;
    %[~,~, target_data_res_thrcov_PCA] = fxn_Marchenko_thrcover_PCA(target_data, 90);
    [Mahalanobis_d_raw{4},  Mahalanobis_mean(:,4)] = fxn_nmt_Mahalanobis(reference_data_PCA_score,target_data_PCA_score.thrcov_PCA_score);% mean of refer data
catch
    warning('wrong');
end
end