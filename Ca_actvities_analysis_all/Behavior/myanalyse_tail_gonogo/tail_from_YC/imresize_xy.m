function xyres = imresize_xy(im, sf, xy)
% xyres = imresize_xy(im, sf, xy)
siz = size(im);
siz = siz(end:-1:1);
ori = (siz + 1) / 2;
orir = (siz*sf+1)/2;
oriall = repmat(ori, size(xy,1), 1);
xyres = (xy - oriall) * sf + repmat(orir, size(xy,1), 1);
% xyres = bsxfun(@plus, bsxfun(@minus, xy, ori) * sf, ori);

