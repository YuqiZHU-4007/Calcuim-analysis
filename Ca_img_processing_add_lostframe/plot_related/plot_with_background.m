function plot_with_background(x, y, bk, plotfunc)
% plot_with_background(x, y, bk, plotfunc)
% plot_with_background(y, bk, plotfunc)
% plot_with_background(y, bk)
% 

if nargin == 3
    if strcmp(class(bk), 'function_handle')
        plotfunc = bk;
        bk = y;
        y = x;
        x = [];
    else
        plotfunc = @plot;
    end
end

if nargin == 2
    plotfunc = @plot;
    bk = y;
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
% 
% figure; 
% plotfunc(x, y);








