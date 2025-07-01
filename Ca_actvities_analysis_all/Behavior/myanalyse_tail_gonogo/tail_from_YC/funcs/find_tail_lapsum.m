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
prolongfac = 2;        % could possibly decrease;

movescope = 30;
movestep  = 3;

% supermask = engim;
x0 = pos0(1); y0 = pos0(2);

posarray = zeros(npoints, 2);
posarray(1,:) = pos0;
thetaarray = zeros(npoints, 1);
thetaarray(1) = theta0;  % upperleft

for jj = 2:npoints
    w = segwidth(jj);
    baseth = thetaarray(jj-1);
    thetas = baseth-movescope:movestep:baseth+movescope;
    
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
        
        siz = ceil(max([r * prolongfac + 9, 2*w+9]));
        [ri, ci] = ind2sub(size(engim), ind);
        tim = zeros(siz);
        tind = sub2ind(size(tim), ri-min(ri)+3, ci-min(ci)+3);
        tim(tind) = engim(ind);
        ls = fspecial('laplacian');
        lsim = imfilter(tim, ls);
        smv = mean(abs(lsim(:)));
%         tim = engim(ind);
        values(ii) = mean(tim(:))  - smv/2;
        
%       figure; imshow(tim); %xlabel(num2str(smv)); ylabel(num2str(values(ii)));

    end
    
    [maxv, maxi] = max(values);
    
    theta = thetas(maxi);
    y = y0 - sind(theta) * r;
    x = x0 + cosd(theta) * r;
    
    posarray(jj,:) = [x y];
    thetaarray(jj) = theta;
    x0 = x; y0 = y;
end



