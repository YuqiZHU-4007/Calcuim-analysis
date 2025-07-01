

p = fileparts(mfilename('fullpath'));

addpath(genpath([p '/../load/']));

addpath(genpath([p '/../utils/']));

addpath(genpath([p '/funcs/']));

%addpath(genpath([p '/plot/']));

%% tools
addpath([p '/..//tools/fiji-win64/Fiji.app/scripts/']);

addpath([p '/../tools/vlfeat/toolbox/']);
vl_setup;

%addpath([p '/../tools/saveastiff/']);
%addpath([p '/../tools/NIfTI/']);
addpath([p '/../tools/icp/']);
addpath([p '/../tools/imagejroi/']);

addpath([p '/../tools/oopsi-master/']);

addpath([p '/../tools/FDA/FDA/']);

addpath(genpath([p '/../tools/utils/']));
addpath(genpath([p '/../tools/formats/']));

addpath(genpath([p '/../tools/randomforest-matlab/']));






