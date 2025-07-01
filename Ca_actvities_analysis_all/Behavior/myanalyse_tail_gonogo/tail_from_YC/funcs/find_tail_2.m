function [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, r, params)

% engim = p(:,:,238);
% engim = imread('med3.tiff', 206);
% [sizey, sizex] = size(engim);
% x0 = 300; y0 = 302;
% x0 = 288; y0 = 297;
% 
% theta0 = 45;  % as right-handed
% npoints = 20;
% r = 15;

npoints = npoints + 1;  % find one more points for smooth

segwidth = linspace(10, 3, npoints);  % tail width  % could possibly decrease;
% w = 7;
prolongfac = 1.5;        % could possibly decrease;

movescope = 30;
movestep  = 5;

% supermask = engim;
x0 = pos0(1); y0 = pos0(2);

posarray = zeros(npoints, 2);
posarray(1,:) = pos0;
thetaarray = zeros(npoints, 1);
thetaarray(1) = theta0;  % upperleft

baseth = thetaarray(1);
thetas = baseth-movescope:movestep:baseth+movescope;

values = zeros(1, length(thetas));
for ii = 1:length(thetas)
    theta = thetas(ii);
    ye = y0 - sind(theta) * (r * prolongfac);
    xe = x0 + cosd(theta) * (r * prolongfac);
    
    ind = roi_line_segment(engim, [x0, y0], [xe ye], segwidth(1));
    values(ii) = mean(engim(ind));
end

for jj = 2:npoints
    baseth = thetaarray(jj-1);
    thetas = baseth-movescope:movestep:baseth+movescope;    
    
%     values = zeros(1, length(thetas));
%     for ii = 1:length(thetas)
%         theta = thetas(ii);
% %         y = y0 - sind(theta) * r;
% %         x = x0 + cosd(theta) * r;
%         ye = y0 - sind(theta) * (r * prolongfac);
%         xe = x0 + cosd(theta) * (r * prolongfac);
%         
%         % insertShape is slow
%         %mask = false([sizey, sizex]);
%         %mask = insertShape(mask, 'line', [x0 y0 x y], 'color', 'w',
%         %'linewidth', w);  
%         %mask = mask(:,:,1) > 0;
%         %mask(ind) = true;
%         %figure; imshow(mask)
%         %mask = imdilate(mask, strel('disk', w));
%         
%         ind = roi_line_segment(engim, [x0, y0], [xe ye], w);
%         values(ii) = mean(engim(ind));
%     end
    
    if jj <= 10
        thres = 0.99;
    else
        thres = 0;
    end
    keepind = find(values >= thres * max(values));
    
    
    values2 = zeros(length(thetas));
    for ii1 = keepind
        theta = thetas(ii1);
        y1 = y0 - sind(theta) * r;
        x1 = x0 + cosd(theta) * r;
        
        thetas2 = theta-movescope:movestep:theta+movescope;
        for ii2 = 1:length(thetas2)
            theta2 = thetas2(ii2);
            y2 = y1 - sind(theta2) * (r * prolongfac);
            x2 = x1 + cosd(theta2) * (r * prolongfac);
            ind2 = roi_line_segment(engim, [x1, y1], [x2 y2], segwidth(jj));
            values2(ii1, ii2) = values(ii1) + mean(engim(ind2));
        end
    end
        
    
    [maxv, maxi] = max(values2(:));
    [ii1, ii2] = ind2sub(size(values2), maxi);
    
    values = values2(ii1,:) - values(ii1);
    
    theta = thetas(ii1);
    y = y0 - sind(theta) * r;
    x = x0 + cosd(theta) * r;
    
    posarray(jj,:) = [x y];
    thetaarray(jj) = theta;
    x0 = x; y0 = y;
end



