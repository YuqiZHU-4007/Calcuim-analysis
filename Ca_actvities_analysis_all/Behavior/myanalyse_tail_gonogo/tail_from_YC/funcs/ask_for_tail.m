function [selxy] = ask_for_tail(frame1, npoints)

%% interactive
% enhframe1 = mat2gray(frame1, vlim);
% enhframe1 = adapthisteq(frame1);
% enhframe1 = frame1;
enhframe1 = normim(frame1, [0, 1/100]);

figure; imshow(enhframe1);
h = impoly();
pos = wait(h);
x = pos(:,1); y = pos(:,2);

if size(pos, 1) >= 3
    windowWidth = 3;
    polynomialOrder = 2;
    smx = sgolayfilt(x, polynomialOrder, windowWidth);
    smy = sgolayfilt(y, polynomialOrder, windowWidth);
else
    smx = x; smy = y;
end
% smx = x; smy = y;
imshow(enhframe1);
hold on;
plot(smx, smy, 'b-*');

%
cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
%taillen = cumlen(end);
selx = interp1([0 cumlen'], smx, linspace(0, cumlen(end), npoints));
sely = interp1([0 cumlen'], smy, linspace(0, cumlen(end), npoints));
% % % seglen = cumlen(end)/(npoints-1);
% % % pos0 = [selx(1) sely(1)];

selxy = [selx(:) sely(:)];

figure; imshow(enhframe1);
hold on;
plot(selx, sely, 'r-*');



