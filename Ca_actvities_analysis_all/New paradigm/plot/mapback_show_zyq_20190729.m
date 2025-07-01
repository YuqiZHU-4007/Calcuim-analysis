function showspv  = mapback_show_zyq_20190729(values, spv, index, mask, drawtype)
% vmap = mapback(values, spv, siz, ind)
%
% values should already be normalized
%
if size(values,3)<3 && ~strcmp(drawtype, 'filledcircle')
    values=repmat(values,1,1,3);
end
if nargin < 4
    index = 1:numel(values);
end
siz=size(mask);
[spvzind, spvregion] = spvindex(spv, siz);

%% map back
if siz(3)>1
    vmap = rescalegd(im2double(mask), [1/10000 1/10000]);
    if strcmp(drawtype, 'dot')
        showspv = zeros(siz(1), siz(2), 3, siz(3));
    elseif strcmp(drawtype, 'circle')
        showspv = zeros(siz(1), siz(2), 3, siz(3));
    end
    for zi = 1:siz(3)
        slice1 = vmap(:,:,zi);slice2 = vmap(:,:,zi);slice3 = vmap(:,:,zi);%slice= zeros(siz(1)*siz(2), 3);
        pti =  index(find(spv(index,3)==zi));
        ind = sub2ind(size(vmap(:,:,zi)), spv(pti, 2), spv(pti, 1));
        if strcmp(drawtype, 'dot')
            for ii=1:length(pti)
            slice1(spvregion{pti(ii)}) = values(find(index==pti(ii)),:,1);  %1;
            slice2(spvregion{pti(ii)}) = values(find(index==pti(ii)),:,2);  %1;
            slice3(spvregion{pti(ii)}) = values(find(index==pti(ii)),3);  %1;
            end
            slice=cat(3,slice1,slice2,slice3);
            %slice = reshape(slice, [siz(1),siz(2), 3]);
            showspv(:,:,:,zi) = slice;
        elseif strcmp(drawtype, 'circle')
            slice = vmap(:,:,zi);
            slice = insertShape(slice, 'circle', [spv(pti, 1), spv(pti, 2), floor(spv(pti, 5)) - 1]);
            showspv(:,:,:,zi) = slice;
        elseif strcmp(drawtype, 'filledcircle')
            slice = vmap(:,:,zi);
            slice = insertShape(slice, 'filledcircle', [spv(pti, 1), spv(pti, 2), floor(spv(pti, 5)) - 1],'color',values(find(spv(index,3)==zi),:));
            showspv(:,:,:,zi) = slice;
        end
    end
end
%mapback



