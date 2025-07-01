function [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, r, params)
% [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, r, params)
%
% find the tail trajectory from tail probability map using greedy opt
%

% %
% engim = p(:,:,238);
% engim = imread('med3.tiff', 206);
% [sizey, sizex] = size(engim);
% x0 = 300; y0 = 302;
% x0 = 288; y0 = 297;
% 
% theta0 = 45;  % as right-handed
% npoints = 20;
% r = 15;
%

if nargin < 6
    params = [];
end
if isfield(params, 'prolong_factor')
    prolongfac = params.prolong_factor;
else
    prolongfac = 2;
end
if isfield(params, 'turn_scope')
    movescope = params.turn_scope;
else
    movescope = 30;
end
if isfield(params, 'turn_step')
    movestep = params.turn_step;
else
    movestep = 2;
end
if isfield(params, 'max_width')
    max_width = params.max_width;
else
    max_width = 11;
end
if isfield(params, 'min_width')
    min_width = params.min_width;
else
    min_width = 3;
end
if isfield(params, 'basetheta')
    basetheta = params.basetheta;
else
    basetheta = NaN;
end

npoints = npoints + 1;  % find one more points for smooth

segwidth = linspace(max_width, min_width, npoints);  % tail width (diameter) % could possibly decrease;

% supermask = engim;
x0 = pos0(1); y0 = pos0(2);

posarray = zeros(npoints, 2);
posarray(1,:) = pos0;
thetaarray = zeros(npoints, 1);
thetaarray(1) = theta0;  % upperleft

for jj = 2:npoints
    w = segwidth(jj);
    baseth = thetaarray(jj-1);
    %thetas = baseth-movescope:movestep:baseth+movescope;
    if jj == 2  % restrain the first point
        th1 = baseth:-movestep:max([basetheta-60, baseth-movescope]);
        th2 = baseth:movestep:min([basetheta+60, baseth+movescope]);
        thetas = [th1(end:-1:1) th2(2:end)];
    else
        thetas = baseth-movescope:movestep:baseth+movescope;
    end
        
    values = zeros(1, length(thetas));
    for ii = 1:length(thetas)
        theta = thetas(ii);
%         y = y0 - sind(theta) * r;
%         x = x0 + cosd(theta) * r;
        ye = y0 - sind(theta) * (r * prolongfac);
        xe = x0 + cosd(theta) * (r * prolongfac);
        
        % insertShape is slow
        %mask = false([sizey, sizex]);
        %mask = insertShape(mask, 'line', [x0 y0 x y], 'color', 'w',
        %'linewidth', w);  
        %mask = mask(:,:,1) > 0;
        %mask(ind) = true;
        %figure; imshow(mask)
        %mask = imdilate(mask, strel('disk', w));
        
        ind = roi_line_segment(engim, [x0, y0], [xe ye], w);
        values(ii) = mean(engim(ind));
    end
    
    [maxv, maxi] = max(values);
    
    theta = thetas(maxi);
    y = y0 - sind(theta) * r;
    x = x0 + cosd(theta) * r;
    
    posarray(jj,:) = [x y];
    thetaarray(jj) = theta;
    x0 = x; y0 = y;
end



