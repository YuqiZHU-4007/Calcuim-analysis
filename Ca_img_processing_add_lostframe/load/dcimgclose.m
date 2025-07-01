function success = dcimgclose(fhandle)
% success = dcimgclose(fhandle)

% % % errmsgp = libpointer('int32Ptr', 0);

success = calllib('tmcamcon', 'TMCC_CLOSEDCIMGFILE', fhandle.fid, fhandle.errmsgp);
if success==0, error('error occured in TMCC_CLOSEDCIMGFILE'); end

