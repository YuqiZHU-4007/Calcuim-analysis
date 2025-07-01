function xlat = frame_register_xlat_prepared_reference_backup_afterward(pre_xlat, mvim)
% xlat = frame_register_xlat_prepared_reference(pre_xlat, targim)

% I1 = double(fxim);
% I1 = I1 - mean(I1(:));
I2 = double(mvim);
I2 = I2 - mean(I2(:));

C = ifft2(pre_xlat .* (fft2(I2)));  % actually not fully understand, but validated as work when movement is small, tm style impossible, 
C = fftshift(C);%将零频点移到频谱的中间
score = max(C(:));
[y, x] = find(C == score);
y = y(1); x = x(1); %% as C is single-layer, just one maxmimum is needed

dx = x - round(size(pre_xlat, 2) / 2 - .5) - 1;%%%
dy = y - round(size(pre_xlat, 1) / 2 - .5) - 1;
xlat = [dx dy];
