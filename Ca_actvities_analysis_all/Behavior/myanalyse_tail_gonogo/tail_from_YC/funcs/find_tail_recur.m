function [posarray, paththeta] = find_tail_recur(engim, pos0, theta0, npoints, r)

% engim = p(:,:,238);
% engim = imread('med3.tiff', 206);
[sizey, sizex] = size(engim);

% theta0 = 45;  % as right-handed
% pos0   = [300, 302];
% pos0   = [290, 297];
% npoints = 15;
% r = 20;

npoints = npoints + 1;  % find two more points for smooth

segwidth = linspace(10, 3, npoints);  % tail width  % could possibly decrease;
% w = 7;
prolongfac = 2;        % could possibly decrease;

movescope = 30;
movestep  = 5;

% w = 7;  % tail width
% movescope = 45;
% movestep  = 10;
prop = 0.8;

[value, paththeta] = recur(theta0, pos0, npoints-1);

% cal posarray from thetaarray
posarray = zeros(npoints, 2);
posarray(1,:) = pos0;
xx0 = pos0(1); yy0 = pos0(2);
for pti = 2:length(paththeta)
    th = paththeta(pti);
    yy = yy0 - sind(th) * r;
    xx = xx0 + cosd(th) * r;
    posarray(pti,:) = [xx yy];
    xx0 = xx;
    yy0 = yy;
end

function [value, theta] = recur(theta0, pos0, level)
    
    x0 = pos0(1);
    y0 = pos0(2);
    
    [candthetas, candvalues] = next(theta0, pos0, segwidth(npoints - level));
    
    if level == 1
        [value, maxi] = max(candvalues);
        theta = candthetas(maxi);
        return;
    end
    
    ncand = length(candthetas);
    pathvalues = zeros(ncand, 1);
    paththetas = zeros(ncand, level);
    for ii = 1:ncand
        theta = candthetas(ii);
        y = y0 - sind(theta) * r;
        x = x0 + cosd(theta) * r;
        [termvalue, termtheta] = recur(theta, [x y], level - 1);
        pathvalues(ii) = candvalues(ii) + termvalue;
        paththetas(ii,:) = termtheta;
    end
    [value, maxi] = max(pathvalues);
    theta = [candthetas(maxi) paththetas(maxi,:)];
end

function [candthetas, candvalues] = next(theta0, pos0, w)

    x0 = pos0(1); y0 = pos0(2);
    thetas = theta0-movescope:movestep:theta0+movescope;
    
    values = zeros(1, length(thetas));
    for ii = 1:length(thetas)
        theta = thetas(ii);
        ye = y0 - sind(theta) * (r * prolongfac);
        xe = x0 + cosd(theta) * (r * prolongfac);
        
        ind = roi_line_segment(engim, [x0, y0], [xe ye], w);
        values(ii) = mean(engim(ind));
    end
    
    maxv = max(values);
    minv = min(values);
    mind = find(values >= minv + prop * (maxv - minv));
    candthetas = thetas(mind);
    candvalues = values(mind);
end

end





