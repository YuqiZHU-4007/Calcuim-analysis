function [frame] = tiffread_zyq_0430(obj, frmid)
% [frame, timestamp] = dcimgread(fhandle, frmid)
%
fhandle=obj.fhandle;
imgpath=  obj.fhandle.imgpath;
% if isempty(obj.t0)
%     obj.t0=0;
% end
ti=floor((frmid)/obj.T)+1;
zi=mod((frmid),obj.T); 
if zi==0 && ti~=0
     zi=obj.T;
end
%frmid=obj.t0 + obj.T * (ti-1) + zi - 1;

frame = imread(fullfile(imgpath,['z',num2str(zi,'%02d')],[num2str(ti,'%04d') '.tif']));
end
%timestamp = fhandle.timestamp.Value;

