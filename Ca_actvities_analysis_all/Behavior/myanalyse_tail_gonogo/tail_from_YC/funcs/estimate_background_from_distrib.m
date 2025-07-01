function [estback, backim] = estimate_background_from_distrib(pixdistrib, vlim, quantx)
% estback = estimate_background_from_distrib(pixdistrib, vlim = [0 255], quantx = 0.05)
%
%
%

% assume the small values form a gaussian distribution and estimate the mean.

% params
if nargin < 3
    quantx = 0.05;
end
if nargin < 2
    vlim = [0 255];
end

% const
npix = size(pixdistrib, 1);
nt   = sum(pixdistrib(1,:));

% normalize
toohigh = sum(pixdistrib(:,vlim(2)+1:end), 2);
pixdistrib(:,vlim(2)+1:end) = 0;
pixdistrib(:,vlim(2)) = pixdistrib(:,vlim(2)) + toohigh; %%keeping total counts constant

% what about too low?

x = 0:255;

e = sum(repmat(x, [npix, 1]) .* pixdistrib/nt, 2);  % mean/expectation of each pixel

% % % figure; imshow(uint8(round(reshape(e, size(frame1)))));

xx = repmat(x, [npix, 1]);
ee = repmat(e, [1, 256]);
pp = pixdistrib/nt;

sd = sqrt(sum((xx - ee).^2 .*  pp, 2));  % sd of each pixel

% % % figure; hist(sd);

noisesd = prctile(sd, 5);  % at least 5% of the pixels are pure background

% quantx = 0.05;  % 
xshift = icdf('norm', quantx, 0, noisesd);  % assume the noise distribution is G(0, noisesd)
if isnan(xshift)
    xshift = 0;
end

% compute quantx quantile
pixcdf = cumsum(pixdistrib, 2) / nt;
backim = zeros(npix, 1);
tic;
for ii = 1:npix
    backim(ii) = find(pixcdf(ii,:) > quantx, 1) - 1;
% % %     incrind = find(diff(pixcdf(ii,:)) > 0);
% % %     %     ls = find(pixcdf(ii,:)==0, 1, 'last');
% % %     %     rs = find(pixcdf(ii,:)==1, 1, 'first');
% % %     if isempty(incrind)
% % %         backim(ii) = pixcdf(ii,1);
% % %     elseif length(incrind) < 2
% % %         backim(ii) = incrind;
% % %     else
% % %         %backim(ii) = interp1(pixcdf(ii,incrind+1), incrind-1, quantx, 'linear', 0);
% % %         backim(ii) = find(pixcdf(ii,:) > quantx) - 1;
% % %     end
end
toc;

% estback = backim;
estback = backim - xshift;

% backim = reshape(backim - xshift, size(frame1));

