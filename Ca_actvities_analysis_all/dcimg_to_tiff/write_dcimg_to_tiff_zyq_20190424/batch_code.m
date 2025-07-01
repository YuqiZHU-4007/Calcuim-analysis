function [pathes]=batch_code(display) 
%batch code -----------------------------
[listfile, dirpath]= uigetfile('G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_nonlearner_cutmove_20190921\path\*.*',display);
if listfile == 0
    displayresult('no file selected');
    error('no file selected');
end
fid = fopen([dirpath '/' listfile]);
pathes = textscan(fid, '%s', 'delimiter', '\n');
pathes = pathes{1};
fclose(fid);