function fhandle = dcimgopen(imgpath)
% fhandle = dcimgopen(imgpath)
%

if ~libisloaded('tmcamcon')
    loadlibrary('tmcamcon.dll', 'tmcamcon.h');
end

% declarations
filep   = libpointer('uint32Ptr', 0);
nframep = libpointer('int32Ptr', 0);
errmsgp = libpointer('int32Ptr', 0);

widthp    = libpointer('int32Ptr', 0);
heightp   = libpointer('int32Ptr', 0);
timestamp = libpointer('doublePtr', 0);
% errmsg  = libpointer('string', '');

success = calllib('tmcamcon', 'TMCC_CLOSEDCIMGFILE', 0, errmsgp);
% no matter whether it will success, try to close.

success = calllib('tmcamcon', 'TMCC_OPENDCIMGFILE', filep, imgpath, nframep, errmsgp);
if success==0, error('error occured in TMCC_OPENDCIMGFILE'); end  
success = calllib('tmcamcon', 'TMCC_GETDCIMGFRAMEINFO', filep.Value, 0, widthp, heightp, errmsgp);
if success==0, error('error occured in TMCC_GETDCIMGFRAMEINFO'); end

framep = libpointer('uint16Ptr', zeros(widthp.Value*heightp.Value, 1));  % in case it's not assigned, easy to find problem

fhandle.fid = filep.Value;
fhandle.nframes = nframep.Value;
fhandle.width   = widthp.Value;
fhandle.height  = heightp.Value;

% cache
fhandle.errmsgp = errmsgp;
fhandle.timestamp = timestamp;
fhandle.framep  = framep;
