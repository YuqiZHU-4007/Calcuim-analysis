function processing_main
global key
global opt
global env
global batchi
global pathes
global pathes_signal
%% average stack
if key==0
    batchi=0;
else
    batchi=batchi+1;
    displayresult(['Starting bactch' num2str(batchi)]);
    sample_path = pathes{batchi};
    displayresult(sample_path);
   
    try
        close(imio);
    catch
    end
    opt=getpara;
    opt.imgpath = sample_path;
    opt.process=struct;
    
    if isempty(pathes_signal)
        signal_path=[];
        displayresult('no signal path');
        opt.process.T_lessframe=[];
        opt.process.z_lessframe=[];
        opt.process.T_moreframe=[];
        opt.process.z_moreframe=[];
    else
        if isempty(pathes_signal{batchi})
            signal_path=[];
            displayresult('no signal path');
            opt.process.T_lessframe=[];
            opt.process.z_lessframe=[];
            opt.process.T_moreframe=[];
            opt.process.z_moreframe=[];
        else
            signal_path=pathes_signal{batchi};
            opt.process.signal_path=signal_path;
            displayresult(signal_path);
            data= xlsread(signal_path,'re_ind_less_signal_use');
            opt.process.T_lessframe=data(2:end,1);
            opt.process.z_lessframe=data(2:end,2:end);
            data= xlsread(signal_path,'re_ind_more_signal_use');
            opt.process.T_moreframe=data(2:end,1);
            opt.process.z_moreframe=data(2:end,2:end);
        end
    end
    
    %[filename, filepath]= uigetfile('F:\fear conditioning\*.xls');
    
    
    %% settings ---------
    [dirpath, filename] = fileparts(opt.imgpath);
    opt.savepath = checkpath([dirpath '/summary_' filename '/' filename '_calcium_raw/']);
    % opt.t0 = 26;
    %% conserved options
    % opt.burstsign = 1e10;  % depricated: maybe reused in burstmode
    depth = max([1, opt.T - 1]);
    %% constants
    % imio = ImageIO(opt.imgpath, opt.T, opt.t0);
    imio = dcimgreader_zyq(opt.imgpath, opt.T, opt.t0,opt.process);  %, opt.t0); % auto estimate t
     displayresult('imio');
    displayresult(imio);
    if imio.nframes == 0  % dll error
        displayresult('There may be error occurred when reading the dcimg file, please restart the Matlab');
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
    displayresult('opt');
    displayresult(opt);
    
    
    %% average volume
    displayresult('computing average volume...');
    disp 'computing average volume...';
    if key==0
    else
        if exist([opt.savepath '/vol.tiff'], 'file') > 0
            warning('reference volume already exists, read from the file!');
            displayresult('warning: reference volume already exists, read from the file!')
            vol = tiffread([opt.savepath '/vol.tiff']);
           % imshow(vol);
        else
            vol = ref_volume_zyq(imio, opt.nsamples, opt.zi,opt.process,opt.savepath);  % caution: inline params
            % use one z slice for motion detection now;%%%%%ÅÜ³öÀ´maskÈ«ºÚ
           % imshow(vol(:,:,1));
            tiffwrite(vol, [opt.savepath '/vol.tiff']);
        end
        env.vol = vol;
    end
    
    %% registration
    displayresult('running registration') ;
    disp 'running registration';
    if key==0
    else
        if opt.isregistration
            if env.depth > 1
                tic; rio = pipeline_registration_ori_zyq(env, opt); 
                displayresult(['run registration time: ' num2str(toc)]);%toc;
            else  % single core registration for single plane
                tic; rio = pipeline_registration_single_thread(env, opt); 
                displayresult(['run registration time: ' num2str(toc)]);%toc;
            end
        else
            %
            temp = 1;
        end  % end if opt.isregistration
    end
    
    %% f-b segmentation         001.0                                         
    displayresult('running background segmentation...') ;
    disp 'running background segmentation...';
    if key==0
    else
        % caution: you may like to check the result of f-b segmentation
        %[env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
        
        if opt.loc_nucleus
            [env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
        else
            [env.volmask, env.mask] = background_segmentation(env, opt);
        end  %%changed 20171026, previously _nucleus used without selection loop
        tiffwrite(im2uint8(env.volmask), [opt.savepath '/volmask.tiff']);
        tiffwrite(im2uint8(cat(2, rescalegd(env.vol, [0 0]), env.volmask)), [opt.savepath '/inspect_volmask.tiff']);
        %volmask = tiffread( [opt.savepath '/volmask.tiff']);
        %env.volmask=volmask;
        %volmask = volmask > 100;
    end
         
    %% supervoxel segmentation
    displayresult('running supervoxel segmentation...') ;
    disp 'running supervoxel segmentation...';
    if key==0
    else
        tic;
        if opt.loc_nucleus
            [env.supervoxel, scoremaps, radmaps] = supervoxel_segmentation_n(env, opt);
        else
            [env.supervoxel, scoremaps, radmaps] = supervoxel_segmentation_cytosol(env, opt);
        end
        toc;
        displayresult(['run supervoxel segmentation time: ' num2str(toc)]);
        %% index supervoxel
       [spvzind, spvregion] = supervoxel_index(env, opt);
%        save([opt.savepath '/supervoxel_index.mat'], 'spvzind', 'spvregion',  '-v7.3');
       save([opt.savepath '/supervoxel_index.mat'], 'spvzind', 'spvregion');
        %% save segmentation result
        showspv = show_spv(env, opt, 'dot');
        svpath = checkpath([opt.savepath '/spvseg_dot/']);
        seqwrite(showspv, svpath);
        showspv = show_spv(env, opt, 'circle');
        svpath = checkpath([opt.savepath '/spvseg_circle/']);
        seqwrite(showspv, svpath);
    end
%     save([opt.savepath '/env.mat'], 'env', 'opt', '-v7.3');
    save([opt.savepath '/env.mat'], 'env', 'opt');
    
    %% ---------------------------------------------------------------
%     %% extract signal
    displayresult('extracting signals...') ;
    disp 'extracting signals....';
    if key==0
    else
        tic; %
        [duration_extract,activities,backgroundy,backPSTH ]=extract_signal(opt,env);
         displayresult(['run extracting signals time: ' num2str(toc)]); toc;  
    end
   % structure2txt(opt,'opt',[opt.process.signal_path(1:end-6) '\opt.txt']);
end

end