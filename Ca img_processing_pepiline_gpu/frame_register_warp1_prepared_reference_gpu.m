function [trans] = frame_register_warp1_prepared_reference_gpu(pre_warp, targvol, warpingSettings)
% [trans] = frame_register_warp_prepared_reference(pre_warp, targim, warpingSettings)
% both inputs are gpuArray, outputing gpuArray

%% registration
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
maximumShift = warpingSettings(3);

[height, width, depth] = size(targvol);

%jVector = (squareSize + 1):pointDensity:(size(image, 1) - squareSize);  % point coord along y
%kVector = (squareSize + 1):pointDensity:(size(image, 2) - squareSize);  % point coord along x
[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));

% target vars

%q = 0;
%for j = jVector
%    for k = kVector
%        q = q + 1;
trans = gpuArray(single(zeros(depth, numel(xx), 2)));

%     for dd = 1:depth
% %         dxArray = gpuArray(zeros(numel(xx),1));  % length(jVector) * length(kVector) = number of control points
% %         dyArray = gpuArray(zeros(numel(xx),1));
%         for q = 1:numel(xx)
%             jj = xx(q); ii = yy(q); 
%             tile2 = targvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), dd); %%already in GPU
%             %         I1 = refrim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
%             tile1fft_gpu = pre_warp(:,:,dd,q); %%sure, fft-ed, confirmed 20170223
%             %I1 = imageSquareReference;
%             %I2 = imageSquare;
%             % todo: caution: why not center the image? 
%             C = ifftshift(ifft2(tile1fft_gpu .* conj(fft2(tile2))));
% %         C = ifftshift(ifft2(tile1fft .* conj(fft2(tile2))));
%         %tic; C = normxcorr2(I1, I2); toc;  % this is way slow
%             [trans(dd,q,2), trans(dd,q,1)] = find(C == max(C(:)));
%             trans(dd,q,1) = trans(dd,q,1) - squareSize - 2;
%             trans(dd,q,2) = trans(dd,q,2) - squareSize - 2;
%             
% %             if length(trans(dd,q,1)) ~= 1
% %                 trans(dd,q,1) = 0;
% %                 trans(dd,q,1) = 0;
% %             end
%             
% %             dxArray(q) = dx;
% %             dyArray(q) = dy;
%         end
%         % todo: now you can add adaptive grid size
%         %dxArray = dxArray .* (abs(dxArray) < maximumShift);
%         %dyArray = dyArray .* (abs(dyArray) < maximumShift);
%         trans(dd,(abs(trans(dd,:,1)) > maximumShift),1) = 0;
%         trans(dd,(abs(trans(dd,:,2)) > maximumShift),2) = 0;
% %         trans (dd, :, :)= [dxArray(:) dyArray(:)];
%     end %%depth
   
%         dxArray = gpuArray(zeros(numel(xx),1));  % length(jVector) * length(kVector) = number of control points
%         dyArray = gpuArray(zeros(numel(xx),1));

blockSize = 10;
Nblock = ceil(numel(xx)/blockSize);
for q_q = 1:Nblock
        blockSize1 = blockSize+(numel(xx)-Nblock*blockSize)*(q_q==Nblock);
        tile2 = gpuArray(single(ones(squareSize*2+1,squareSize*2+1,depth,blockSize1)));        
        for q = 1:blockSize1
            jj = xx(q+(q_q-1)*blockSize); ii = yy(q+(q_q-1)*blockSize); 
            tile2(:,:,:,q) = targvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), :); %%already in GPU
            %         I1 = refrim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
%             tile1fft_gpu = pre_warp(:,:,:,q); %%sure, fft-ed, confirmed 20170223
            %I1 = imageSquareReference;
            %I2 = imageSquare;
            % todo: caution: why not center the image? 

%          for q = 1:numel(xx)
%             for dd = 1:depth               
% %                 [trans(dd,q,2), trans(dd,q,1)] = find(C(:,:,dd,q) == C_max(:,:,dd,q));
%                 trans(dd,q,1) = (trans(dd,q,1) - squareSize - 2)*((trans(dd,q,1) - squareSize - 2)<=maximumShift);
%                 trans(dd,q,2) = (trans(dd,q,2) - squareSize - 2)*((trans(dd,q,2) - squareSize - 2)<=maximumShift);
%             end
%          end
            
%             if length(trans(dd,q,1)) ~= 1
%                 trans(dd,q,1) = 0;
%                 trans(dd,q,1) = 0;
%             end
            
%             dxArray(q) = dx;
%             dyArray(q) = dy;
        end
        C_mult = pre_warp(:,:,:,(q_q-1)*blockSize+(1:blockSize1)) .* conj(fft2(tile2));
        C = ifftshift(ifft2(C_mult));
        clear tile2 C_mult
        C_max = repmat(max(max(C)),squareSize*2+1,squareSize*2+1,1,1);
%       C = ifftshift(ifft2(tile1fft .* conj(fft2(tile2))));
        %tic; C = normxcorr2(I1, I2); toc;  % this is way slow
        C_ind = find(C == C_max);
        [C_ind2,C_ind1,C_ind3,~] = ind2sub(size(C),C_ind);
        C_ind3_diff = [nan;diff(C_ind3)];
        C_ind2 = C_ind2(C_ind3_diff~=0) -squareSize -2;
        C_ind1 = C_ind1(C_ind3_diff~=0) -squareSize -2;
        trans1 = reshape([C_ind1 C_ind2],depth,blockSize1,2); %corrected20170417
        trans1 (abs(trans1)>maximumShift) =0; %%corrected 20170417
        clear C C_max C_ind2 C_ind1 C_ind
        trans(:,(q_q-1)*blockSize+(1:blockSize1), :) = trans1;
        clear trans1
        % todo: now you can add adaptive grid size
        %dxArray = dxArray .* (abs(dxArray) < maximumShift);
        %dyArray = dyArray .* (abs(dyArray) < maximumShift);
%         trans(dd,(abs(trans(dd,:,1)) > maximumShift),1) = 0;
%         trans(dd,(abs(trans(dd,:,2)) > maximumShift),2) = 0;
%         trans (dd, :, :)= [dxArray(:) dyArray(:)];
end

end

