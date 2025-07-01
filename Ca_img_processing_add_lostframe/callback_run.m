function callback_run
clear all;
global key
global opt
global env
global batchi
global pathes
global pathes_signal

key=1;
hfc=findobj('tag','hmainfigure');
hdis=findobj(hfc,'tag','dispanel');
set(findobj(hdis,'tag','dis_tag'),'string',[]);
key=1;
if ~isempty( timerfind('tag','testtimer'))%testtimer�Ѿ������У�ֹͣ��ɾ����ǰtesttimer
    stop(timerfind('tag','testtimer'));
    delete(timerfind('tag','testtimer'));
end
%batch code -----------------------------
[listfile, dirpath]= uigetfile('z:/*.*');
if listfile == 0
    displayresult('no file selected');
    error('no file selected');
end
fid = fopen([dirpath '/' listfile]);
pathes = textscan(fid, '%s', 'delimiter', '\n');
pathes = pathes{1};
nbatch = length(pathes);batchi=0;
fclose(fid);

[signalfilename, signalfilepath]= uigetfile('z:/*.*');
if signalfilename == 0
    displayresult('no file selected');
    warning('no file selected');
    pathes_signal=[];
else
    fid = fopen([signalfilepath '/' signalfilename]);
pathes_signal = textscan(fid, '%s', 'delimiter', '\n');
pathes_signal = pathes_signal{1};
fclose(fid);
end

%% setup dir
file_masterdir =  'F:\DUlab\FC_analyse';%'Z:/';
code_dir = fullfile(file_masterdir,'Ca_img_processing_add_lostframe');
if ~exist(code_dir,'dir')
    error(['can not find' code_dir]);
end
addpath(genpath(code_dir));
cd(code_dir)
setenv('MW_MINGW64_LOC','F:\DUlab\FC_analyse\Ca_img_processing_add_lostframe\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64');

testtimer=timer('executionmode','fixedSpacing','period',1,'TasksToExecute', nbatch,'tag','testtimer');
set(testtimer,'timerfcn',@(x,y)processing_main);
set(testtimer,'startfcn',{@timer_fcn,'start timer'});
set(testtimer,'stopfcn',{@timer_fcn,'stop timer'});
start(testtimer);

end% batch for
% end batch code --




