function xlat = frame_register_xlat_prepared_reference_gpu(pre_xlat, mvvol)
% xlat = frame_register_xlat_prepared_reference(pre_xlat, targim)

% I1 = double(fxim);
% I1 = I1 - mean(I1(:));
% I2 = double(mvvol);
I2 = mvvol - mean(mvvol(:));

% I2_gpu=gpuArray(I2);
C_gpu = ifft2(pre_xlat .* (fft2(I2))); % actually not fully understand, but validated as work when movement is small, tm style impossible
clear I2 
C_gpu = fftshift(C_gpu);
% C = gather(C_gpu);
xlat = gpuArray(single(zeros(size(mvvol,3),1,2)));
        
C_max = repmat(max(max(C_gpu)),size(C_gpu,1),size(C_gpu,2),1);
C_ind = find (C_gpu(:) == C_max(:));
[C_ind2, C_ind1, C_ind3] = ind2sub(size(C_gpu),C_ind);
C_ind3_diff = [nan;diff(C_ind3)];
C_ind2 = C_ind2(C_ind3_diff~=0) -round(size(pre_xlat, 2) / 2 - .5) - 1;
C_ind1 = C_ind1(C_ind3_diff~=0) -round(size(pre_xlat, 1) / 2 - .5) - 1;

xlat(:,1,1) = C_ind1;
xlat(:,1,2) = C_ind2;
% 
% for dd=1: size(mvvol,3)
% score = max(max(C_gpu(:,:,dd)));
% [y, x] = find(C_gpu(:,:,dd) == score);
% 
% xlat(dd,1) = x - round(size(pre_xlat, 2) / 2 - .5) - 1;
% xlat(dd,2) = y - round(size(pre_xlat, 1) / 2 - .5) - 1;
% end

% xlat = [dx dy];
end