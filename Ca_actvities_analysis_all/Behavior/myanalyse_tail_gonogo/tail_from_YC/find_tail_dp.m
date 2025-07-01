function tail = find_tail_dp(engim, xy0, taillen)

% suppose image are rotated
% 
% taillen as hamminglen

x0 = round(xy0(1));
y0 = round(xy0(2));

mx = size(engim,2) - x0;
my = size(engim,1) - y0;
[xx, yy] = meshgrid(0:mx, 0:my);

%% distmap
disttoori = xx + yy;
[sdist, order] = sort(disttoori(:));
nsee = find(sdist > taillen, 1);
xxx = xx(order);
yyy = yy(order);

scoremap = zeros(size(disttoori));
scoremap(1, 1) = engim(y0, x0);
pathmap = zeros(size(disttoori));  % 1 means comes horizontal, 2 means comes vertical

for ii = 2:nsee
    xi = xxx(ii);
    yi = yyy(ii);
    
    if xi == 0
        scorel = -Inf;
    else
        scorel = scoremap(yi+1, xi);
    end
    if yi == 0
        scoreu = -Inf;
    else
        scoreu = scoremap(yi, xi+1);
    end
    
    [s,p] = max([scorel scoreu]);  % 1 for x, 2 for y
    pathmap(yi+1, xi+1) = p;
    scoremap(yi+1, xi+1) = engim(yi+y0, xi+x0) + s;
end

ind = find(disttoori == taillen);
[mv, mi] = max(scoremap(ind));
tind = ind(mi);
[ty, tx] = ind2sub(size(scoremap), tind);
ty = ty - 1; tx = tx - 1;

% figure; imshow(normim(frame))
% figure; imshow(normim(engim));
% hold on; plot(tx, ty, '*r')

traceback = zeros(taillen, 2);
traceback(1,:) = [tx, ty];

k = 2;
xi = tx; yi = ty;
while xi ~= 0 || yi ~= 0
    p = pathmap(yi+1, xi+1);
    
    if p == 1
        xi = xi - 1;
    elseif p == 2
        yi = yi - 1;
    else
        error('unexpected');
    end
    traceback(k,:) = [xi yi];
    k = k + 1;
end

tail = traceback;
tail(:,1) = traceback(:,1) + x0;
tail(:,2) = traceback(:,2) + y0;
tail = tail(end:-1:1,:);
% figure; imshow(normim(engim));
% hold on; plot(traceback(:,1), traceback(:,2), 'r')
% 
% hold on; plot(selx, sely, 'y');

