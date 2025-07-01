function ind = roi_line_segment(im, p1, p2, linewidth)
% ind = roi_line_segment(im, p1, p2, linewidth)
% linewidth refers to radius
%

siz = size(im);
pos = round([p1; p2]);

lb = max([1 round(min(pos(:,1)) - linewidth)]);
rb = min([siz(2) round(max(pos(:,1)) + linewidth)]);
tb = max([1 round(min(pos(:,2)) - linewidth)]);
bb = min([siz(1) round(max(pos(:,2)) + linewidth)]);
[xx, yy] = meshgrid(lb:rb, tb:bb);
xx = xx(:); yy = yy(:);
pt = [xx yy];

base = sqrt(sum(abs(p1 - p2).^2));
side1 = dist(pt, p1');
side2 = dist(pt, p2');

semiperi = (base + side1 + side2) / 2;

s = sqrt(semiperi .* (semiperi - base) .* (semiperi - side1) .* (semiperi - side2));
d = 2 * s ./ base; %%distance from one point to the line connecting two others points, called Heron formula

segind = find(d <= linewidth & base*base >= abs(side1.^2 - side2.^2));
ri = yy(segind);
ci = xx(segind);

ind = sub2ind(siz, ri, ci);


