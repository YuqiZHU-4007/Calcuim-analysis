function traceback = find_tail_dp(engim, selxy)


selx = selxy(:,1);  sely = selxy(:,2);

x0 = round(selx(1));
y0 = round(sely(1));

cumlen = cumsum(sqrt(diff(selx).^2 + diff(sely).^2));
taillen = cumlen(end);
nv = [selx(2) - selx(1), sely(2) - sely(1)];


[xx, yy] = meshgrid(1:size(engim,2), 1:size(engim,1));
% distmap = abs(xx - x0) + abs(yy - y0);

%% valid region
xx = xx - x0; yy = yy - y0;
aa = nv(1) * xx + nv(2) * yy < 0;
bb = xx.^2 + yy.^2 > taillen*taillen;
xx = xx + x0; yy = yy + y0;

% figure; imagesc(aa);
% figure; imagesc(bb);

temp = ones(size(engim)-2);
temp = padarray(temp, [1 1]) > 0;
validregion = ~aa&~bb & temp;
% figure; imagesc(validregion)

%% distmap
disttoori = abs(xx(:) - x0) + abs(yy(:) - y0);
ind = 1:length(disttoori);

ind = ind(validregion);
disttoori = disttoori(validregion);

[disttoori, order] = sort(disttoori);
ind = ind(order);
xxx = xx(ind);
yyy = yy(ind);

scoremap = zeros(size(engim));
scoremap(y0, x0) = engim(y0, x0);
pathmap = zeros(size(engim));  % 1 means comes horizontal, 2 means comes vertical

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

hamminglen = floor(abs(sely(end)- sely(1)) + abs(selx(end) - selx(1)));
% hamminglen = ceil(taillen);

% [mv, mi] = max(scoremap(distmap==hamminglen));

ii = find(disttoori == hamminglen);
[mv, mi] = max(scoremap(ind(ii)));
tind = ind(ii(mi));
[ty, tx] = ind2sub(size(engim), tind);

% figure; imshow(normim(frame))
% figure; imshow(normim(engim));
% hold on; plot(tx, ty, '*r')

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
    k = k + 1;
end

% figure; imshow(normim(engim));
% hold on; plot(traceback(:,1), traceback(:,2), 'r')
% 
% hold on; plot(selx, sely, 'y');

