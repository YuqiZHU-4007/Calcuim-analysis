clear all;
load('F:\5.fear conditioning behavioral data\data\20210709\fish2\activities_aft_process.mat')
load('F:\5.fear conditioning behavioral data\data\20210709\fish2\behavior_from_Results_of_alltheta.mat')

% Cal, ��ÿһ�н����ز���
A = activities_preCS_dfdf_aftcorrect; % ʾ������
factor = 2.5; % ����������������
p = 2; % 
q = 5; % 
Cal_downsampled = zeros(size(A, 1)/factor, size(A, 2));
for i = 1:size(A, 2)
   Cal_downsampled(:, i) = resample(A(:, i), p, q);
end
Cal_downsampled=Cal_downsampled';
% % % ȥ�쳣ֵ
% % = Cal_downsampled;
% % Q1 = prctile(data, 25); % ���ķ�λ��
% % Q3 = prctile(data, 75); % ���ķ�λ��
% % IQR = Q3 - Q1; % �ķ�λ��
% % lowerBound = Q1 - 1.5 * IQR; % �½�
% % upperBound = Q3 + 1.5 * IQR; % �Ͻ�
% % dataCleaned = data;
% % dataCleaned(data < lowerBound | data > upperBound) = NaN; % ���쳣ֵ��Ϊ NaN
% % 
% % [xx,yy]=find(data>1);
% % positions = [xx, yy];
% % disp(positions);

% BehQav
x = y_3sd;
M = 60; % ����������
Behav_downsampled = downsample(x, M); % ����������ź�

%% plot
baseColormap = [
  0 0 0;       % ��ɫ
    1 0 0;       % ��ɫ
    1 1 0];      % ��ɫ
numPoints = 5000;

% ��ֵ��չ��ɫӳ��
customColormap = interp1(linspace(1, 3, 3), baseColormap, linspace(1, 3, numPoints), 'linear');

figure;
imagesc(Cal_downsampled);
colormap(customColormap);
colorbar
caxis([0 0.2]);

figure;
plot(Behav_downsampled,'color','k');
xlim([0 1080]);






