function [freq,freqpeaknum,freqmax,freqmaxloc]=findfreq(test)
l=length(test);
freq=abs(fft(test))*2/l;%normalization
freq=freq(1:floor(l/2));
if length(freq)<3
    locs=1;
else
    [pks,locs]=findpeaks(freq);
end
freqpeaknum=length(locs);
% if length(pks)==0
%     [m,loc]=max(freq);
%     freqmax=m;
%     freqmaxloc=loc;
% else
%     [m,loc]=max(pks);
%     freqmax=m;
%     freqmaxloc=locs(loc);
% end
[freqmax,freqmaxloc]=max(freq);
end