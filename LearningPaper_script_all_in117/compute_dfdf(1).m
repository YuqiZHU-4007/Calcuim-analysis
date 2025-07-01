function [dfdf, baselines] = compute_dfdf(ts, winsize, prct)
% dfdf = compute_dfdf(ts, winsize = 33, prct = 0.2)
% compute $\delta F / F_0$ using percentile filter
% 
% inputs:
%   ts: calcium time series
%   winsize: window size used in ordfilt2
%   5*Fs okay for photo-evoked response
%

if nargin < 3
    prct = 0.2;
end

%boffset = 0.1 * mean(ts(:));
if nargin < 2
    winsize = 33;
end

if mod(winsize, 2) == 0
    error('winsize must be odd!');
end

prop = prct;
hwin = (winsize - 1) / 2;
order = floor(prop * winsize);
baselines = ordfilt2(ts, order, true(1, winsize));
boffset = mean(baselines(:)) / 10;

dfdf = (ts - baselines) ./ (baselines + boffset);
%dfdf = (ts - baselines) ./ (baselines);
%dfdf(dfdf<0) = 0;
dfdf = dfdf(:,hwin+1:end-hwin);

tail = repmat(mean(dfdf, 2), [1, hwin]);
dfdf = [tail dfdf tail];
% tail = repmat(mean(ba, 2), [1, hwin]);
baselines(:,1:hwin) = repmat(baselines(:,hwin+1), 1, hwin);
baselines(:,end-hwin+1:end) = repmat(baselines(:,end-hwin), 1, hwin);
end

%{
env.ts;
baselines = dfdf;
for svi = 1:size(env.ts, 1)
    baseline = smooth(env.ts(svi,:), 21)';  % caution: inline parameters
    baselines(svi,:) = baseline;
    dfdf(svi,:) = (env.ts(svi,:) - baseline) ./ baseline;
end
dfdf(isnan(dfdf)) = 0;
%env.dfdf = dfdf;
%}
