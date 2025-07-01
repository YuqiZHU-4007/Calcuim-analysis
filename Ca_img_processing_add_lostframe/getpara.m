%%%%%%%%%%%%%%%%%%%%%把当前面板参数调入para全局结构体%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%已改，str2double忘记%%%%%%%%%%%%%%%%%%%%%%
function curpara=getpara

hfc=findobj('tag','hmainfigure');
hopt=findobj(hfc,'tag','optpanel');
hdis=findobj(hfc,'tag','dispanel');
%hdisedit=findobj(hdis,'tag','dis_tag');

heditT=findobj(hopt,'tag','T_tag');%1
heditt0=findobj(hopt,'tag','t0_tag');%2
% heditT2=findobj(hopt,'tag','T2_tag');%3
heditpoolsize=findobj(hopt,'tag','poolsize_tag');%2

heditloc_nucleus=findobj(hopt,'tag','loc_nucleus_tag');%4
heditisregistration=findobj(hopt,'tag','isregistration_tag');%5
heditextractsignal=findobj(hopt,'tag','extractsignal_tag');%6
heditnsample=findobj(hopt,'tag','nsample_tag');%7
heditzi=findobj(hopt,'tag','zi_tag');%8
heditminrad=findobj(hopt,'tag','minrad_tag');%9
heditmaxrad=findobj(hopt,'tag','maxrad_tag'); %10
heditthres=findobj(hopt,'tag','thres_tag');%11
heditdensity=findobj(hopt,'tag','density_tag');%12
heditsquaresize=findobj(hopt,'tag','squaresize_tag');%13
heditshift=findobj(hopt,'tag','shift_tag');%14


curpara.T=str2double(get(heditT,'string'));
curpara.t0=str2double(get(heditt0,'string'));
curpara.poolsize=str2double(get(heditpoolsize,'string'));
%curpara.T2=str2double(get( heditT2,'string'));

curpara.loc_nucleus=(get(heditloc_nucleus,'value'));
curpara.isregistration=(get(heditisregistration,'value'));
curpara.extract_signal_multicore=(get(heditextractsignal,'value'));
curpara.nsamples=str2double(get(heditnsample,'string'));
curpara.zi=str2double(get(heditzi,'string'));
curpara.minrad=str2double(get(heditminrad,'string'));
curpara.maxrad=str2double(get(heditmaxrad,'string'));
curpara.thres=str2double(get(heditthres,'string'));

curpara.minrad=str2double(get(heditminrad,'string'));
curpara.maxrad=str2double(get(heditmaxrad,'string'));
curpara.density=str2double(get(heditdensity,'string'));
curpara.squaresize=str2double(get(heditsquaresize,'string'));
curpara.shift=str2double(get(heditshift,'string'));

curpara.warpopt=[curpara.density curpara.squaresize curpara.shift];
if isnan(curpara.t0) curpara.t0=[];end
if isnan(curpara.thres) curpara.thres=[];end
if curpara.loc_nucleus==2 
    curpara.loc_nucleus=false;
else
    curpara.loc_nucleus=true;
end
if curpara.isregistration==2 
    curpara.isregistration=false;
else
    curpara.isregistration=true;
end
if curpara.extract_signal_multicore==2
    curpara.extract_signal_multicore=false;
else
    curpara.extract_signal_multicore=true;
end
end