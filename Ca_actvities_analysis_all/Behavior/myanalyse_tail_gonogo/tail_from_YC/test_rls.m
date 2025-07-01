
% use fitting for find tails
% failed

selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);

% initial values.
theta0 = atan2d(sely(2)-sely(1), selx(2)-selx(1)) - atan2d(0, 1);  % [-180, 180]
theta0 = - theta0;  % find_tail use right-handed coord
thetap = theta0;
pos0   = [selx(1) sely(1)];
cumlen = cumsum(sqrt(diff(selx).^2 + diff(sely).^2));
taillen = cumlen(end);
seglen = cumlen(end)/(npoints-1);
params.basetheta = theta0;

x0 = round(selx(1));
y0 = round(sely(1));
% % % normalvec = [selx(2) - selx(1), sely(2) - sely(1)];
% % % x = 0:0.1:1;
% % % %normalvec(1) * x + normalvec(2) * y = 0
% % % y = normalvec(1) * x / normalvec(2);
% % % 
% % % figure; plot()

[xx, yy] = meshgrid(1:size(frame1,2), 1:size(frame1,1));
xx = xx - x0;
yy = yy - y0;

x = 0:1:100;
nv = [selx(2) - selx(1), sely(2) - sely(1)];
y = - normalvec(1) * x / normalvec(2);
figure; imshow(normim(frame1))
hold on; plot(selx - x0, sely - y0, 'r'); 
hold on; plot(x, y, 'g')
ylim([-500 500])

xx.^2 + yy.^2 > taillen || nv(1) * xx + nv(2) * yy < 0
aa = nv(1) * xx + nv(2) * yy < 0;
bb = xx.^2 + yy.^2 > taillen*taillen;

frame(aa|bb) = 0;
engim(aa|bb) = 0;
% 
% % temp = engim / sum(engim(:)) * 100000;
% % t2 = round(temp);
% % t2(t2==1) = 0;
% % 
% % pts = image2density(t2);
% % pts(:,1) = pts(:,1) - x0;
% % pts(:,2) = pts(:,2) - y0;
% % 
% % figure; plot(pts(:,1), pts(:,2), '.')
% % 
% % ft = fittype('a*x + b*x^2 + c*x^3', 'coeff', {'a', 'b', 'c'});
% % obj = fit(pts(:,1), pts(:,2), ft, 'robust', 'Bisquare');
% % figure; plot(obj, pts(:,1), pts(:,2))

frame = x023329;
frame = double(frame);
frame = frame ./ sum(frame(:));
frame = normim(frame, [0, 1/100]);
frame = round(frame .* 255);

engim = abs(frame - backim);

selx = selxy(:,1);  sely = selxy(:,2);


resize = 11;
engimr = imresize(engim, 1/resize, 'box');
selxr = (selx - size(frame1,2)/2) / resize + size(engimr,2)/2;
selyr = (sely - size(frame1,1)/2) / resize + size(engimr,1)/2;

figure; imshow(normim(engimr));
hold on; plot(selxr, selyr, 'r');
selx = selxr;
sely = selyr;
engim = engimr;

x0 = round(selx(1));
y0 = round(sely(1));
[xx, yy] = meshgrid(1:size(engim,2), 1:size(engim,1));
distmap = abs(xx - x0) + abs(yy - y0);
disttoori = abs(xx(:) - x0) + abs(yy(:) - y0);
ind = 1:length(disttoori);

cumlen = cumsum(sqrt(diff(selx).^2 + diff(sely).^2));
taillen = cumlen(end);
nv = [selx(2) - selx(1), sely(2) - sely(1)];
xx = xx - x0; yy = yy - y0;
aa = nv(1) * xx + nv(2) * yy < 0;
bb = xx.^2 + yy.^2 > taillen*taillen;
xx = xx + x0; yy = yy + y0;

figure; imagesc(aa);
figure; imagesc(bb);

temp = ones(size(engim)-2);
temp = padarray(temp, [1 1]) > 0;
validregion = ~aa&~bb & temp;
figure; imagesc(validregion)

ind = ind(validregion);
disttoori = disttoori(validregion);

[disttoori, order] = sort(disttoori);
ind = ind(order);
xxx = xx(ind);
yyy = yy(ind);

scoremap = zeros(size(engim));
scoremap(y0, x0) = engim(y0, x0);
pathmap = zeros(size(engim));  % 1 means comes horizontal, 2 means comes vertical

tic;
for ii = 2:length(ind)
    xi = xxx(ii);
    yi = yyy(ii);
    
    temp = xi - x0;
    xpu = temp / abs(temp);
    temp = yi - y0;
    ypu = temp / abs(temp);
    if isnan(xpu), xpu = 0; end
    if isnan(ypu), ypu = 0; end
    
    [s,p] = max([scoremap(yi, xi-xpu), scoremap(yi-ypu, xi)]);
    pathmap(yi, xi) = p;
    scoremap(yi, xi) = engim(yi, xi) + s;
end
toc;

% hamminglen = round(abs(sely(end)- sely(1)) + abs(selx(end) - selx(1)));
hamminglen = ceil(taillen);

% [mv, mi] = max(scoremap(distmap==hamminglen));

ii = find(disttoori == hamminglen);
[mv, mi] = max(scoremap(ind(ii)));
tind = ind(ii(mi));
[ty, tx] = ind2sub(size(engim), tind);

% figure; imshow(normim(frame))
figure; imshow(normim(engim));
hold on; plot(tx, ty, '*r')

traceback = zeros(hamminglen, 2);
traceback(1,:) = [tx, ty];

k = 2;
xi = tx; yi = ty;
while xi ~= x0 || yi ~= y0
    p = pathmap(yi, xi);
    
    % prev unit
    temp = xi - selx(1);
    xpu = temp / abs(temp);
    temp = yi - sely(1);
    ypu = temp / abs(temp);
    if isnan(xpu), xpu = 0; end
    if isnan(ypu), ypu = 0; end
    
    if p == 1
        xi = xi - xpu;
    elseif p == 2
        yi = yi - ypu;
    else
        error('');
    end
    traceback(k,:) = [xi yi];
    [xi yi]
    k = k + 1;
end

figure; imshow(normim(engim));
hold on; plot(traceback(:,1), traceback(:,2), 'r')

hold on; plot(selx, sely, 'y');






















