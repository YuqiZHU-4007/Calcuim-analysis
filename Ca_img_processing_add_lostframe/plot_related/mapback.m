function vmap = mapback(values, spv, siz, ind)
% vmap = mapback(values, spv, siz, ind)
%
% values should already be normalized
%

if nargin < 4
    ind = 1:numel(values);
end

[spvzind, spvregion] = spvindex(spv, siz);

%% map back
if siz(3)>1
    vmap = zeros(siz(1)*siz(2), siz(3));
    for ii = 1:numel(ind)
        spi = ind(ii);
        vmap(spvregion{spi}, spv(spi,3)) = values(ii);
    end
    vmap = reshape(vmap,siz);%max(vmap(:))
else
    vmap = zeros(siz(1)*siz(2), 3);
    for ii = 1:numel(ind)
        spi = ind(ii);
        vmap(spvregion{spi}, 1) = values(ii,1);
        vmap(spvregion{spi}, 2) = values(ii,2);
        vmap(spvregion{spi}, 3) = values(ii,3);
    end
    vmap = reshape(vmap, [siz(1),siz(2), 3]);
end

% seqwrite(vmap, svpath);


