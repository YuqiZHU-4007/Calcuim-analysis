function scaledMatrix = scaleMatrixToRange(matrix, newMin, newMax)
% scaleMatrixToRange �������ֵ���ŵ�ָ���ķ�Χ
% matrix: �������
% newMin: ���ź����Сֵ
% newMax: ���ź�����ֵ

% �ҵ���������ֵ����Сֵ
maxVal=max(matrix(:));
minVal =min(matrix(:));

% �����������Ӻ�ƫ����
scaleFactor = (newMax - newMin) / (maxVal - minVal);
offset = newMin - scaleFactor * minVal;

% Ӧ�����ź�ƫ��
scaledMatrix = scaleFactor * matrix + offset * ones(size(matrix));
end
