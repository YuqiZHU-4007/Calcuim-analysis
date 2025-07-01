function pts = image2density(im)

npt = sum(im(:));
pts = size(npt, 2);

[y, x] = find(im);

pti = 1;
for ii = 1:length(x)
    ni = im(y(ii),x(ii));
    pts(pti:pti+ni-1,1) = x(ii);
    pts(pti:pti+ni-1,2) = y(ii);
    pti = pti + ni;
end

pts = pts + randn(size(pts))/2;


