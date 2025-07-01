function [pathes]=batch_code(display) 
%batch code -----------------------------
[listfile, dirpath]= uigetfile('C:\standard brain\*.*',display);
if listfile == 0
    displayresult('no file selected');
    error('no file selected');
end
fid = fopen([dirpath '/' listfile]);
pathes = textscan(fid, '%s', 'delimiter', '\n');
pathes = pathes{1};
fclose(fid);