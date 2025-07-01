function [moving] = estimate_threshold_for_movement_detection(stdwin, winsize, basethres)

if nargin < 3
    basethres = 1e-10;
end

stdwin = abs(stdwin);

hwin = ceil((winsize-1)/2);
st = hwin + 1;
ed = find(stdwin(st:end) < basethres, 1);
if isempty(ed)
    ed = length(stdwin);
else
    ed = ed - hwin;
end
s = stdwin(st:ed);

bs = env_secant(1:length(s), s, winsize * 30, 'bottom');

figure; plot(s); hold on; plot(bs, 'r');

ss = s(:) - bs(:);
ss = ss';

figure; plot(ss);

candv = linspace(min(ss), max(ss), 100);

areas = zeros(size(candv));
for ii = 1:length(candv)
    thres = candv(ii);
    areas(ii) = sum(ss > thres);
end

%figure; plot(candv, areas)

areav = diff(areas);

minthres = prctile(ss, 25); 

[~, locs] = findpeaks(areav);

ii = find(candv(locs) > minthres, 1);

thres = candv(locs(ii));
if isempty(thres) 
    warning('estimate threshold failed!');
    thres = minthres * 2;
end

figure; plot(candv(2:end), areav); 
hold on; plot(candv(locs(ii)), areav(locs(ii)), '^r');

figure; plot(ss); 
hold on; plot(1:length(ss), repmat(thres, size(ss)), 'r');

moving = ss > thres;
lab = bwlabel(moving);
counts = hist(lab, unique(lab));

%figure; bar(counts(2:end));

moving(counts < hwin) = false;
temp = false(size(stdwin));
temp(st:ed) = moving;
moving = temp;

%movind = find(moving);
%movind = movind + st - 1;

figure; plot(stdwin); hold on; plot(moving * max(stdwin), 'r');


