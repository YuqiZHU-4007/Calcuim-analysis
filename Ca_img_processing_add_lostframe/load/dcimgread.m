function [frame, timestamp] = dcimgread(fhandle, frmid)
% [frame, timestamp] = dcimgread(fhandle, frmid)
%

% % % timestamp = libpointer('doublePtr', 0);
% % % errmsgp = libpointer('int32Ptr', 0);
% % % framep = libpointer('uint16Ptr', zeros(fhandle.width*fhandle.height, 1));
% pre-alloc have no speed difference, but theoritically, pre-alloc.
% if isnan(frmid)
%     frame = nan(fhandle.height, fhandle.width);
% else
%success = 
calllib('tmcamcon', 'TMCC_GETDCIMGFRAMEDATA_A', fhandle.fid, frmid-1, fhandle.framep, fhandle.timestamp, fhandle.errmsgp);
% % % calllib('tmcamcon', 'TMCC_GETDCIMGFRAMEDATA_A', fhandle.fid, frmid-1, framep, timestamp, errmsgp);
% if success==0, error(['error occured in TMCC_GETDCIMGFRAMEDATA_A when reading the ' num2str(frmid) 'th frame.']); end
frame = reshape(fhandle.framep.Value, fhandle.width, fhandle.height)';
end
%timestamp = fhandle.timestamp.Value;

