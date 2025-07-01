function [spvzind, spvregion] = spv_index(env, opt)
% [spvzind, spvregion] = spv_index(env, opt)

spv = env.supervoxel;
height = env.height;
width = env.width;
depth = env.depth;

spvzind = cell(depth,1);
for zi = 1:depth
    spvzind{zi} = find(spv(:,3) == zi)';
end

nspv = size(spv, 1);
spvregion = cell(nspv, 1);
for pi = 1:nspv
    xp = spv(pi,1);
    yp = spv(pi,2);
    rp = floor(spv(pi,4));  % use outer radius?
    [x, y] = meshgrid(xp-rp-1:xp+rp+1,yp-rp-1:yp+rp+1);
    ind = find((x - xp).^2 + (y - yp).^2 <= rp.^2);
    r = y(ind); c = x(ind);
    try
        ind = sub2ind([height width], r, c);%sub2ind函数将矩阵中指定元素的行列下标转换成存储的序号，即线性索引号。
    catch
        r(r<1) = 1; r(r>height) = height;
        c(c<1) = 1; c(c>width) = width;
        ind = sub2ind([height width], r, c);
        ind = unique(ind);
    end
    spvregion{pi} = ind;
end

