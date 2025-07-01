function [posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen, params)
% [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, seglen, params)
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

npoints = npoints + 1;  % find one more points for smooth

%totallen = npoints * seglen;
segwidth = linspace(max_width, min_width, npoints);  % tail width (radius) % could possibly decrease;
%seglen   = totallen / sum(linspace(2, 1, npoints)) .* linspace(2, 1, npoints);

% supermask = engim;
% object
posarray = zeros(npoints, 2);
posarray(1,:) = pos0;
thetaarray = zeros(npoints, 1);
thetaarray(1) = thetap;  % upperleft

%% find point 2 and 3 at first.
x0 = pos0(1); y0 = pos0(2);
w = segwidth(2);
r = seglen;
thetas1 = thetap-movescope:movestep:thetap+movescope;
values = zeros(length(thetas1));
for ii = 1:length(thetas1)
    theta1 = thetas1(ii);
    ye = y0 - sind(theta1) * (r * prolongfac);
    xe = x0 + cosd(theta1) * (r * prolongfac);
    y1 = y0 - sind(theta1) * r;
    x1 = x0 + cosd(theta1) * r;
    
    ind = roi_line_segment(engim, [x0, y0], [xe ye], w);
    values(ii,:) = mean(engim(ind));
    
    w = segwidth(3);
    %r = seglen;
    thetas2 = theta1-movescope:movestep:theta1+movescope;
    for jj = 1:length(thetas2)
        theta2 = thetas2(jj);
        y2 = y1 - sind(theta2) * (r * prolongfac);
        x2 = x1 + cosd(theta2) * (r * prolongfac);
        ind = roi_line_segment(engim, [x1, y1], [x2 y2], w);
        values(ii,jj) = values(ii,jj) + mean(engim(ind));
    end 
end

[~, ind] = max(values(:));
[ii, jj] = ind2sub(size(values), ind);
theta1 = thetas1(ii);
y1 = y0 - sind(theta1) * r;
x1 = x0 + cosd(theta1) * r;
thetas2 = theta1-movescope:movestep:theta1+movescope;
theta2 = thetas2(jj);
y2 = y1 - sind(theta2) * r;
x2 = x1 + cosd(theta2) * r;
posarray(2:3,:) = [x1 y1; x2 y2];
thetaarray(2:3) = [theta1; theta2];

x0 = x2; y0 = y2;
%% find remaining tails
for jj = 4:npoints
    w = segwidth(jj);
    %r = seglen(jj);
    
    baseth = thetaarray(jj-1);
    thetas = baseth-movescope:movestep:baseth+movescope;
    
    values = zeros(1, length(thetas));
    for ii = 1:length(thetas)
        theta = thetas(ii);
%         y = y0 - sind(theta) * r;
%         x = x0 + cosd(theta) * r;
        ye = y0 - sind(theta) * (r * prolongfac);
        xe = x0 + cosd(theta) * (r * prolongfac);
                
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



