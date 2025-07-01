clear all;
load('F:\5.fear conditioning behavioral data\data\20210709\fish2\activities_aft_process.mat')
load('F:\5.fear conditioning behavioral data\data\20210709\fish2\behavior_from_Results_of_alltheta.mat')

% Cal, 对每一列进行重采样
A = activities_preCS_dfdf_aftcorrect; % 示例矩阵
factor = 2.5; % 非整数降采样因子
p = 2; % 
q = 5; % 
Cal_downsampled = zeros(size(A, 1)/factor, size(A, 2));
for i = 1:size(A, 2)
   Cal_downsampled(:, i) = resample(A(:, i), p, q);
end
Cal_downsampled=Cal_downsampled';
% % % 去异常值
% % = Cal_downsampled;
% % Q1 = prctile(data, 25); % 下四分位数
% % Q3 = prctile(data, 75); % 上四分位数
% % IQR = Q3 - Q1; % 四分位距
% % lowerBound = Q1 - 1.5 * IQR; % 下界
% % upperBound = Q3 + 1.5 * IQR; % 上界
% % dataCleaned = data;
% % dataCleaned(data < lowerBound | data > upperBound) = NaN; % 将异常值设为 NaN
% % 
% % [xx,yy]=find(data>1);
% % positions = [xx, yy];
% % disp(positions);

% BehQav
x = y_3sd;
M = 60; % 降采样因子
Behav_downsampled = downsample(x, M); % 降采样后的信号

%% plot
baseColormap = [
  0 0 0;       % 黑色
    1 0 0;       % 红色
    1 1 0];      % 黄色
numPoints = 5000;

% 插值扩展颜色映射
customColormap = interp1(linspace(1, 3, 3), baseColormap, linspace(1, 3, numPoints), 'linear');

figure;
imagesc(Cal_downsampled);
colormap(customColormap);
colorbar
caxis([0 0.2]);

figure;
plot(Behav_downsampled,'color','k');
xlim([0 1080]);






