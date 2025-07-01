function [shake_individual_mean_std] = shake_individual_bin(shake_all,trial)
%����ÿ������ѧϰǰ�󶯵����
%   shake_individual_mean_std:������ά�ȵĵ�һά�Ǿ�ֵ���ڶ�ά�Ǳ�׼��

shake_individual_mean_std = [];
shake_individual_mean_std(1,:,1) = mean(shake_all(1:trial.hab,:),1,'omitnan');
shake_individual_mean_std(1,:,2) = std(shake_all(1:trial.hab,:),1,1,'omitnan');
for ii = 1:trial.block
    shake_individual_mean_std(ii+1,:,1) = mean(shake_all((trial.hab+trial.perblock*(ii-1)+1):(trial.hab+trial.perblock*ii),:),1,'omitnan');
    shake_individual_mean_std(ii+1,:,2) = std(shake_all((trial.hab+trial.perblock*(ii-1)+1):(trial.hab+trial.perblock*ii),:),1,1,'omitnan');
end
shake_individual_mean_std(7,:,1) = mean(shake_all((trial.total-trial.test+1):trial.total,:),1,'omitnan');
shake_individual_mean_std(7,:,2) = std(shake_all((trial.total-trial.test+1):trial.total,:),1,1,'omitnan');

end

