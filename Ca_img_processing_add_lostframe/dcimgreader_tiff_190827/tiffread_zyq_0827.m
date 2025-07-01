function [frame] = tiffread_zyq_0827(obj, frmid)
% [frame, timestamp] = dcimgread(fhandle, frmid)
%20190827,�Ӷ�ȡ����tif
fhandle=obj.fhandle;
imgpath=  obj.fhandle.imgpath;
type=obj.fhandle.filetypeF;
% if isempty(obj.t0)
%     obj.t0=0;
% end
ti=floor((frmid)/obj.T)+1;
zi=mod((frmid),obj.T);
if zi==0 && ti~=0
    zi=obj.T;
end
[a,b]=fileparts(imgpath);
%frmid=obj.t0 + obj.T * (ti-1) + zi - 1;
switch type
    case 'single file'
        if~exist(fullfile(imgpath,[num2str(frmid,'%04d') '.tif']))
            frame = imread(fullfile(imgpath,[b,'_',num2str(frmid,'%01d') '.tif']));
        else
            frame = imread(fullfile(imgpath,[num2str(frmid,'%04d') '.tif']));
        end
    case 'multi file'
        frame = imread(fullfile(imgpath,['z',num2str(zi,'%02d')],[num2str(ti,'%04d') '.tif']));
end
end
%timestamp = fhandle.timestamp.Value;

