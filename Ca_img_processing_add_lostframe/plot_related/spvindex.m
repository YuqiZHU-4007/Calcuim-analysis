function [spvzind, spvregion] = spvindex(spv, siz)
% [spvzind, spvregion] = spvindex(spv, siz)
%

if nargin < 2
    xmax = max(spv(:,1));
    ymax = max(spv(:,2));
    zmax = max(spv(:,3));
    rmax = max(spv(:,5));
    height = ymax + rmax;
    width  = xmax + rmax;
    depth  = zmax;
    
else
    height = siz(1);
    width = siz(2);
    depth = siz(3);

end

spvzind = cell(depth,1);
for zi = 1:depth
    spvzind{zi} = find(spv(:,3) == zi)';
end

nspv = size(spv, 1);
spvregion = cell(nspv, 1);
for pi = 1:nspv
    xp = spv(pi,1);
    yp = spv(pi,2);
    rp = floor(spv(pi,5));  % use outer radius?
    [x, y] = meshgrid(xp-rp-1:xp+rp+1,yp-rp-1:yp+rp+1);
    ind = find((x - xp).^2 + (y - yp).^2 <= rp.^2);
    r = y(ind); c = x(ind);
    try
        ind = sub2ind([height width], r, c);
    catch
        r(r<1) = 1; r(r>height) = height;
        c(c<1) = 1; c(c>width) = width;
        ind = sub2ind([height width], r, c);
        ind = unique(ind);
    end
    spvregion{pi} = ind;
end

