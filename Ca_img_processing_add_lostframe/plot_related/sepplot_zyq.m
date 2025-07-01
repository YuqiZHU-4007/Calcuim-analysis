function [h, incre,sig] = sepplot_zyq(x, y,textind,prct)
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

if nargin == 1
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
incre = repmat((0:interval:interval*nsig-interval)', [1, nt]);
sig = sig + repmat((0:interval:interval*nsig-interval)', [1, nt]);
h = plot(x, sig','b');
str=textind{1};
text(min(x)*ones(size(incre(:,1))),incre(:,1),str,'HorizontalAlignment','right');
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w','Position',[0.0395,-0.0329,0.9484,1.0482],...
    'OuterPosition',[-0.119568158168574,-0.174375739963405,1.22372528616025,1.286190937466364]) ;
%ylabel(str)














