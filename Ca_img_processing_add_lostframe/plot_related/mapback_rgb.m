function vmap = mapback(values, spv, siz, ind)
% vmap = mapback(values, spv, siz, ind)
% 
% values should already be normalized
% 

if nargin < 4
    ind = 1:size(values,1);
end

[spvzind, spvregion] = spvindex(spv, siz);


vset = unique(values)
if any(vset < 0)
    error('expecting index image');
end
[~, values] = histc(values, vset);
nval = length(vset);
if any(vset == 0)
    values = values - 1;
    nval = nval - 1;
end
values;

%% map back
vmap = zeros(siz(1)*siz(2), siz(3));
for ii = 1:size(values,1)
    spi = ind(ii);
    vmap(spvregion{spi}, spv(spi,3)) = values(ii);
end

vmap = ind2rgb(vmap + 1, [0 0 0; hsv(nval)]);

%vmap = reshape(vmap, [siz(1) siz(2) siz(3) 3]);



