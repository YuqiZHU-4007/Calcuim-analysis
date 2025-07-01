function [h, incre] = sepplot(x, y, c,prct)
% [h, baselines] = sepplot(x, y, prct = 95)
% x and y must be vectors, not scalar
%

if nargin == 3
    if isscalar(y)
        prct = y;
        y = x;
        x = [];
    else
        prct = 99.99;
    end
end

if nargin == 2
    prct = 99.99;
    y = x;
    x = [];
end

if isempty(x)
    if isvector(y)
        nx = length(y);
    else
        nx = size(y, 1);
    end
    x = 1:nx;
end

sig = y';

[nsig, nt] = size(sig);

interval = prctile(sig(:), prct);
if interval==0
    interval=10;
    warning('interval is wrong');
end
incre = -repmat((0:interval:interval*nsig-interval)', [1, nt]);
sig = sig - repmat((0:interval:interval*nsig-interval)', [1, nt]);
h = plot(x, sig','LineWidth',1);
for ii=1:length(h)
    h(ii).Color=c(ii,:);
end
ylim([min(sig(:)) max(sig(:))]);
%figure,plot(y)











