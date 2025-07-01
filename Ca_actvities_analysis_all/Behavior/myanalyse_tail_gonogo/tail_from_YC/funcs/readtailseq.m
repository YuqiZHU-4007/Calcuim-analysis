function [frame, nframes] = readtailseq(seqpath, frameid)
% frame = readtailseq(seqpath, frameid)
% frameid starts from 1
% 

fid = fopen(seqpath);
nframes = fread(fid, 1, 'ulong');
headers = fread(fid, 9, 'int');
width = headers(3);
height = headers(4);
framesize = headers(1);
% if nargout == 2
fileinfo = dir(seqpath);
nframes = fileinfo.bytes / framesize - 1;  % never trust the header about the number of frames
% end

% % % 
% % % 
% % % if nframes == 0
% % %     warning('nframe in header is wrong, using file size!');
% % %     fileinfo = dir(seqpath);
% % %     nframes = fileinfo.bytes / framesize - 1;
% % % %     if floor(nframes) ~= nframes
% % % %         error('cannot compute nframes from file size!');
% % % %     end
% % % end
% % %
% % %

if frameid > nframes
    error('frameid exceeds number of frames!');
end
    
fseek(fid, 0, -1);  % header occupy the whole first frame
fseek(fid, framesize * frameid, 0);

pix = fread(fid, framesize, 'uint8=>uint8');
pix = pix(1:width*height);
frame = reshape(pix, width, height);
frame = frame';

fclose(fid);

