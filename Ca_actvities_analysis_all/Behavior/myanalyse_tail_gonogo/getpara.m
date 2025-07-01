function [para]=getpara
%Maximum tail curvature
%
global t_dur_cell
global alltailpos_everytrial
global taildeg_everytrial
global trailnum
k=1;
para=struct();
para.curv={};para.curv1={};
para.length=12000;
for i=1:trailnum
    tdur=t_dur_cell{i,1};
    pos=alltailpos_everytrial{i,1};
    %test=zeros;
    for j=1:size(tdur,1)
        if tdur(j,2)<=para.length
            para.length=tdur(j,2);
        end
        para.marker(k,:)=[i,j];
        para.ang10raw(k,1:tdur(j,2))=taildeg_everytrial(i,tdur(j,1):tdur(j,3));
        para.ang10mean(k,:)=mean(abs(taildeg_everytrial(i,tdur(j,1):tdur(j,3))));
        para.ang10max(k,:)=max(abs(taildeg_everytrial(i,tdur(j,1):tdur(j,3))));%'Maximum tail angle'
        
        %%w/p
        if tdur(j,2)<=2
             [angpeak,angpeak_loc]=max(taildeg_everytrial(i,tdur(j,1):tdur(j,3)));
        else
        [angpeak,angpeak_loc]=findpeaks(taildeg_everytrial(i,tdur(j,1):tdur(j,3)));
        end
        para.angpeak(k,1:length(angpeak_loc))=angpeak;
        para.angpeakloc(k,1:length(angpeak_loc))=angpeak_loc;
        para.angpeaknum(k,:)=length(angpeak_loc);
        a=find(diff(sign(diff(taildeg_everytrial(i,tdur(j,1):tdur(j,3)))))>0)+1;
        b=find(diff(sign(diff(taildeg_everytrial(i,tdur(j,1):tdur(j,3)))))<0)+1;
        IndMin(k,1:length(a))=taildeg_everytrial(i,a+tdur(j,1));%获得局部最小值的位置??
        IndMax(k,1:length(b))=taildeg_everytrial(i,b+tdur(j,1));%获得局部最大值的位置??
        para.ang10peakmaxr(k,:)=max(IndMin(k,:));
        para.ang10peakmaxl(k,:)=min(IndMax(k,:));
        para.ang10r(k,:)=length(find(IndMin(k,1:length(a))<0));
        para.ang10l(k,:)=length(find(IndMax(k,1:length(b))>0));
        %%%
        [para.freq(k,1:floor(tdur(j,2)/2)),para.freqpeaknum(k,:),para.freqmax(k,:),para.freqmaxloc(k,:)]=findfreq(para.ang10raw(k,1:tdur(j,2)));%频率最大点和值
        %%%
        para.tailpos(:,:,1:tdur(j,2))=pos(:,:,tdur(j,1):tdur(j,3));
        %%%curvature
        [para.curv{k},para.curv1{k}]=findcurv(para.tailpos(:,:,1:tdur(j,2)));
        [para.curvmax(k,1:tdur(j,2)),para.curvmaxloc(k,1:tdur(j,2))]=max(abs(para.curv{k}));%tip deflection
        [para.curvmaxmax(k,:),l]=max(para.curvmax(k,1:tdur(j,2)));
        para.curvmaxmaxloc(k,:)=para.curvmaxloc(k,l);
        
        [para.curvmax1(k,1:tdur(j,2)),para.curvmaxloc1(k,1:tdur(j,2))]=max(abs(para.curv1{k}));%tip deflection
        [para.curvmaxmax1(k,:),l]=max(para.curvmax1(k,1:tdur(j,2)));
        para.curvmaxmaxloc1(k,:)=para.curvmaxloc1(k,l);
        %a=tabulate(para.curv_max_loc(k,1:tdur(j,2)));[m,l]=max(a(:,2));
        %para.curv_max_max_loc(k,:)=a(l,1);
        para.dur(k,:)=tdur(j,2);
        k=k+1;
        %plot((0:floor(a(j,2)/2)-1)/fs,fourier(1:floor(a(j,2)/2)));hold on;
        %scatter(locs(loc),para.frenquency(i,j));hold on;
    end
end
end

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

function [curv,curv1]=findcurv(rawdata)
%插值or不插值
% n=200;
% for i=1:size(rawdata,3)
% rawdata(1:n,1,i)=interp1(rawdata(:,2,i),rawdata(:,1,i),linspace(min(rawdata(:,1,i)),max(rawdata(:,1,i)),n),'pchip');
% rawdata(1:n,2,i)=interp1(rawdata(:,1,i),rawdata(:,2,i),linspace(min(rawdata(:,2,i)),max(rawdata(:,2,i)),n),'pchip');
% end
%pointi-pionti+1  i+1curv
vec=rawdata(1:end-1,:,:)-rawdata(2:end,:,:);
curv(:,:)=atan2d(vec(:,2,:),vec(:,1,:));
curv1(:,:)=atand(vec(:,2,:)./vec(:,1,:));
end
