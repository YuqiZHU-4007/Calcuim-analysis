function [para]=getpara_event_20190521(x,ref_win,win,fs)

dim=size(x);

base=mean(reshape(x(ref_win,:),[],1));
smoothwin=3;

for ii=1:dim(2)%col
    xx=x(win,ii);
    xxs=smoothdata(xx,'movmean',smoothwin);
    [m,indm]=max(abs(xxs));
    para.Max_amp(ii)=xxs(indm); 
    para.Max_amp_Point(ii)=indm/fs;%max amp
    half_maximum=(max(abs(xxs))-base)/2 ;     
    ind=find(abs(xxs)>=abs(half_maximum));ind_rise=ind(1);ind_decay=ind(end);

    para.rise_HM(ii)=[indm-ind_rise]/fs;para.rise_HM_Point(ii)=ind_rise/fs; %FWHM 上升相起始点横坐标值
    para.decay_HM(ii)=[ind_decay-indm]/fs;para.decay_HM_Point(ii)=ind_decay/fs; %FWHM 下降相起始点横坐标值
    para.FWHM(ii)=[ind_decay-ind_rise]/fs;  % FWHM full width at half maximum
    para.AUC(ii)=trapz(xxs(ind_rise:ind_decay)); %area
    
end

