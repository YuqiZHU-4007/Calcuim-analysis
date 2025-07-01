function scaledMatrix = scaleMatrixToRange(matrix, newMin, newMax)
% scaleMatrixToRange 将矩阵的值缩放到指定的范围
% matrix: 输入矩阵
% newMin: 缩放后的最小值
% newMax: 缩放后的最大值

% 找到矩阵的最大值和最小值
maxVal=max(matrix(:));
minVal =min(matrix(:));

% 计算缩放因子和偏移量
scaleFactor = (newMax - newMin) / (maxVal - minVal);
offset = newMin - scaleFactor * minVal;

% 应用缩放和偏移
scaledMatrix = scaleFactor * matrix + offset * ones(size(matrix));
end
