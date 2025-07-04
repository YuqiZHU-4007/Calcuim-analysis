%%demo script
%refine from F:\DUlab\FC analyse\Ca_img_processing\processing_main(dcimgreader to dcimgreader_tiff)
%%zyq,20190610
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!! please do the follow things before you run this script:
% change the file_masterdir which is the code folder path
% change listfile and signalfilename which are the txt file contain path of
% processing data and corresponding lost frame data(usually is empty)

% change the opt.T every time according to your number of z-slices.
% change the opt.loc_nucleus everytime, true for nucleus localization and false for cytosol localization
% optional:
%   change the value of opt.isregistration and opt.extract_signal_multicore
%   if you want to change the neuron segmentation parameters: minrad, maxrad and thres
%       for further information, please refer to function supervoxel_segmentation_n
%   if you want to change the registration parameters: opt.warpopt
%       for further information, please refer to function frame_registration_prepared_reference
%   if you want to change the opt.t0 which is auto estimated in default
% notes:
% auto estimate t0 in default, always get rid of the first volume.
% registration is parallel for z slices, single core for single plane.
% extra notes:
% check the caution part
% do not to fully trust every pipeline function
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup dir
file_masterdir = 'F:\DUlab\FC analyse\';% 'F:\DUlab\FC analyse\';
code_dir = fullfile(file_masterdir,'Ca_img_processing_tiff_brief');
if ~exist(code_dir,'dir')
    error(['can not find' code_dir]);
end
addpath(genpath(code_dir));
cd(code_dir)
addpath(genpath('../dcimgreader_tiff'));

% batch code -----------------------------
%batch code -----------------------------
listfile='F:\DUlab\FC analyse\Ca_img_processing_tiff_brief\DATA/data.txt';
fid = fopen(listfile);
pathes = textscan(fid, '%s', 'delimiter', '\n');
pathes = pathes{1};
nbatch = length(pathes);
fclose(fid);

signalfilename='F:\DUlab\FC analyse\Ca_img_processing_tiff_brief\DATA/dataframe.txt';
fid = fopen(signalfilename); 
pathes_signal = textscan(fid, '%s', 'delimiter', '\n');
pathes_signal = pathes_signal{1};
fclose(fid);

opt.poolsize=12;
opt.T = 25;
opt.t0 = 1;  % [] for auto estimate t0, and will be 1 for single plane data.
opt.loc_nucleus = true;
opt.isregistration = true;
opt.extract_signal_multicore = true;

% segmentation and registration params
opt.nsamples =5;  % sample size used for compute the averaged volume
opt.zi=7; % sample plane used for compute the averaged volume
opt.minrad = 3; % minimal radius of neurons
opt.maxrad = 7; % maximal radius of neurons
opt.thres = [];  % [] for estimate %1/2048/2;  % for nucleus
opt.warpopt = [80 128 12]; %density/square size/threshold

for batchi = 1:nbatch
    %% settings ---------
    sample_path = pathes{batchi};
    display(sample_path);
    try
        close(imio);
    catch
    end
    opt.imgpath = sample_path;
    [dirpath, filename] = fileparts(opt.imgpath);
    opt.savepath = checkpath([dirpath '/summary_' filename '/' filename '_calcium_raw/']);
    opt.process=struct;
    if isempty(pathes_signal)
        signal_path=[];
        display('no signal path');
        opt.process.T_lessframe=[];
        opt.process.z_lessframe=[];
        opt.process.T_moreframe=[];
        opt.process.z_moreframe=[];
    else
        if isempty(pathes_signal{batchi})
            signal_path=[];
            display('no signal path');
            opt.process.T_lessframe=[];
            opt.process.z_lessframe=[];
            opt.process.T_moreframe=[];
            opt.process.z_moreframe=[];
        else
            signal_path=pathes_signal{batchi};
            opt.process.signal_path=signal_path;
            display(signal_path);
            data= xlsread(signal_path,'re_ind_less_signal_use');
            opt.process.T_lessframe=data(2:end,1);
            opt.process.z_lessframe=data(2:end,2:end);
            data= xlsread(signal_path,'re_ind_more_signal_use');
            opt.process.T_moreframe=data(2:end,1);
            opt.process.z_moreframe=data(2:end,2:end);
        end
    end
    depth = max([1, opt.T - 1]);
    %% constants
    % imio = ImageIO(opt.imgpath, opt.T, opt.t0);
    imio = dcimgreader_zyq_tiff(opt.imgpath, opt.T, opt.t0,opt.process);  %, opt.t0); % auto estimate t
    if imio.nframes == 0  % dll error
        error('There may be error occurred when reading the dcimg file, please restart the Matlab');
    end
    opt.t0 = imio.t0;
    
    frame = readframe(imio, 1, 1);
    
    [height, width] = size(frame);
    nvol = imio.nt;
    env.imio = imio;
    env.height = height;
    env.width = width;
    env.depth = depth;
    env.nvol = nvol;
    %% debug options
    opt.subsamplez = 1:env.depth;
    opt.subsamplet = 1:env.nvol;
    opt
    
    %% QC: plot sth to see
    % zi = 28;
    % sumframe = zeros(1, env.nvol);
    % for ti = 1:env.nvol
    %     ti
    %     frame = readframe(env.imio, ti, zi);
    %     sumframe(ti) = sum(double(frame(:)));
    % end
    
    %%%%%%%%  average volume
    disp 'computing average volume...';
    if exist([opt.savepath '/vol.tiff'], 'file') > 0
        warning('reference volume already exists, read from the file!');
        vol = tiffread([opt.savepath '/vol.tiff']);
        % imshow(vol);
    else
        vol = ref_volume_zyq(imio, opt.nsamples, opt.zi,opt.process,opt.savepath);  % caution: inline params
        % use one z slice for motion detection now;
        % imshow(vol(:,:,1));
        tiffwrite(vol, [opt.savepath '/vol.tiff']);
    end
    env.vol = vol;
    
    %% registration
    disp 'running registration';
    if opt.isregistration
        if env.depth > 1
            tic; rio = pipeline_registration_ori_zyq_tiff(env, opt);
            toc;
        else  % single core registration for single plane
            tic; %rio = pipeline_registration_single_thread(env, opt);
            error('not contain pipeline_registration_single_thread');
            toc;
        end
    else
        %
        temp = 1;
    end  % end if opt.isregistration
    
    %% f-b segmentation
    disp 'running background segmentation...';
    % caution: you may like to check the result of f-b segmentation
    %[env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
    if opt.loc_nucleus
        [env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
    else
        [env.volmask, env.mask] = background_segmentation(env, opt);
    end  %%changed 20171026, previously _nucleus used without selection loop
    tiffwrite(im2uint8(env.volmask), [opt.savepath '/volmask.tiff']);
    tiffwrite(im2uint8(cat(2, rescalegd(env.vol, [0 0]), env.volmask)), [opt.savepath '/inspect_volmask.tiff']);
    
    %% supervoxel segmentation
    disp 'running supervoxel segmentation...';
    tic;
    if opt.loc_nucleus
        [env.supervoxel, scoremaps, radmaps] = supervoxel_segmentation_n(env, opt);
    else
        [env.supervoxel, scoremaps, radmaps] = supervoxel_segmentation_cytosol(env, opt);
    end
    toc;
    %% index supervoxel
    [spvzind, spvregion] = supervoxel_index(env, opt);
    save([opt.savepath '/supervoxel_index.mat'], 'spvzind', 'spvregion',  '-v7.3');
    %% save segmentation result
    showspv = show_spv(env, opt, 'dot');
    svpath = checkpath([opt.savepath '/spvseg_dot/']);
    seqwrite(showspv, svpath);
    showspv = show_spv(env, opt, 'circle');
    svpath = checkpath([opt.savepath '/spvseg_circle/']);
    seqwrite(showspv, svpath);
    save([opt.savepath '/env.mat'], 'env', 'opt', '-v7.3');
    
    %% extract signal
    disp 'extracting signals....';
    tic; %
    [duration_extract,activities,backgroundy,backPSTH ]=extract_signal_tiff(opt,env);
    toc;
    % batch code ------
end % batch for
% end batch code --



