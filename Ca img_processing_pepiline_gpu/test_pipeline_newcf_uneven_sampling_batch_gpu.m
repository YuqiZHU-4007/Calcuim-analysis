

% batch code -----------------------------
[listfile, dirpath]= uigetfile('z:/*.*');
if listfile == 0
    error('no file selected');
end

nsteps_og_vol = [];
ind_beg_step=[2 2 2 8 8 8 8 8 8 7 4 4 4 6 6 6 6 6 3 3 3 3 6 4 5 8 0 4 5 8 0];
ind_end_step=[4 4 4 10 10 10 10 10 10 9 6 6 6 8 8 8 8 8 5 5 5 5 8 6 7 10 2 6 7 10 2];

fid = fopen([dirpath '/' listfile]);
pathes = textscan(fid, '%s', 'delimiter', '\n');
pathes = pathes{1};
nbatch = length(pathes);
fclose(fid);



for batchi = [5 6 4]
    sample_path = pathes{batchi}
    % end batch code -------------------------
    
    try
        close(imio);
    catch
    end
    
    opt.imgpath = sample_path;
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!!! please do the follow things before you run this script:
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % change the opt.T every time according to your number of z-slices.
    % change the opt.loc_nucleus everytime, true for nucleus localization and false for cytosol localization
    %
    % optional:
    %   change the value of opt.isregistration and opt.extract_signal_multicore
    %
    %   if you want to change the neuron segmentation parameters: minrad, maxrad and thres
    %       for further information, please refer to function supervoxel_segmentation_n
    %   if you want to change the registration parameters: opt.warpopt
    %       for further information, please refer to function frame_registration_prepared_reference
    %   if you want to change the opt.t0 which is auto estimated in default
    %
    
    % notes:
    % auto estimate t0 in default, always get rid of the first volume.
    % registration is parallel for z slices, single core for single plane.
    %
    
    % extra notes:
    % check the caution part
    % do not to fully trust every pipeline function
    
    %% settings ---------
    %opt.imgpath = 'Z:\data\20160420\fish 2/rec122119095.dcimg';
    [dirpath, filename] = fileparts(opt.imgpath);
    filename
    opt.savepath = checkpath([dirpath '/summary_' filename '/' filename '_calcium_raw/']);
    % opt.imgpath = [dirpath '/' imgname];
    % opt.imgpath = 'Z:\data\20160420\fish 3/rec126155512.dcimg';
    % opt.imgpath = 'Z:\data\20160420\fish 1/rec104720631.dcimg';
    % opt.imgpath = 'Z:\data\20160419\fish2/rec25083993.dcimg';
    % opt.imgpath = 'Z:\data\20160419\fish1/rec12013948.dcimg';
    % opt.imgpath = 'Z:\data\20160415/rec8498965.dcimg';
    % opt.t0 = 26;
    opt.T = 14;
    opt.t0 = 1;  % [] for auto estimate t0, and will be 1 for single plane data and 1 for non-frame-losing data.
    
    opt.nvols_per_trial =235;
    opt.ind_og_vol = 160; %% counted from 0 in labview, 1 smaller than real index
    opt.nsteps_stand_vol = opt.T;
    opt.nsteps_og_vol = 72;
    opt.ind_beg_step = ind_beg_step(batchi)+1; %% counted from 0 in labview, 1 smaller than real index
    opt.ind_end_step = ind_end_step(batchi)+1; %% counted from 0 in labview, 1 smaller than real index
    opt.nvolets_og = opt.nsteps_og_vol/(opt.ind_end_step-opt.ind_beg_step)/2;
    opt.nvolets_per_trial = opt.nvols_per_trial-1+opt.nvolets_og;
    opt.nsteps_per_trial = opt.nsteps_stand_vol*(opt.nvols_per_trial-1) + opt.nsteps_og_vol;
    
    % todo: compute t0 from data.
    opt.loc_nucleus = true;
    opt.isregistration = true;     % true
    opt.extract_signal_multicore = true;
    
    
    % segmentation params
    opt.minrad = 5;
    opt.maxrad = 13;
    % opt.thres = 2/256;
    opt.thres = [];  % [] for estimate %1/2048/2;  % for nucleus
    % 50 ~ 70 k, prior known as 100 k
    % 100 k not ga
    %% ------setting over
    
    % savepath
    %opt.savepath = checkpath('Z:/data/20151111/rec72098716_sig/');
    % [dirpath, filename] = fileparts(opt.imgpath);
    % opt.savepath = checkpath([dirpath '/summary_' filename '/' filename '_calcium_raw/']);
    
    % % % opt.imgpath = 'Z:\data\20151110\rec9858264.dcimg';
    % % % opt.t0 = 25;
    % % % opt.T = 25;
    % % % opt.savepath = checkpath('Z:\data\20151110/rec9858264_sig/');
    % % %
    % % % opt.imgpath = 'Z:\data\20151111\testdcimg.dcimg';
    % % % opt.t0 = 25;
    % % % opt.T = 25;
    % % % opt.savepath = checkpath('Z:\data\20151110/testdcimg_sig/');
    
    %% conserved options
    % opt.burstsign = 1e10;  % depricated: maybe reused in burstmode
    depth = max([1, opt.T - 1]);
    
    % registration params
    opt.warpopt = [80 128 12];
    %opt.warpopt = [112 128 10];
    
    % compute t0
    
    %% constants
    % imio = ImageIO(opt.imgpath, opt.T, opt.t0);
    imio = dcimgreader_uneven(opt.imgpath, opt);  %, opt.t0); % auto estimate t0
    
    imio
    
    if imio.nframes == 0  % dll error
        error('There may be error occurred when reading the dcimg file, please restart the Matlab');
    end
    opt.t0 = imio.t0;
    
    frame = readframe(imio, opt, 1, 1);
    
    [height, width] = size(frame);
    nvolet = imio.nt;
    
    env.imio = imio;
    env.height = height;
    env.width = width;
    env.depth = depth;
    env.nvolet = nvolet;
    
    %% debug options
    opt.nsamples = 128;  % sample size used for compute the averaged volume
    opt.subsamplez = 1:env.depth;
    opt.nTrial = floor(double(imio.nframes - opt.t0 + 1) / opt.nsteps_per_trial);
    ind_stand_vol = [1:opt.ind_og_vol opt.ind_og_vol+opt.nvolets_og+1:opt.nvolets_per_trial];
    opt.ind_stand_volet_total = repmat(ind_stand_vol(:),1,opt.nTrial) + repmat(0:opt.nvolets_per_trial:(opt.nTrial-1)*opt.nvolets_per_trial,opt.nvols_per_trial-1,1);
    opt.ind_stand_volet_total = opt.ind_stand_volet_total(:);    
    opt.subsamplet = opt.ind_stand_volet_total;
    
    opt
    
    %% QC: plot sth to see
    % zi = 28;
    % sumframe = zeros(1, env.nvol);
    % for ti = 1:env.nvol
    %     ti
    %     frame = readframe(env.imio, ti, zi);
    %     sumframe(ti) = sum(double(frame(:)));
    % end
    
    %% average stack
    disp 'computing average volume...';
    if exist([opt.savepath '/vol.tiff'], 'file') > 0
        warning('reference volume already exists, read from the file!');
        vol = tiffread([opt.savepath '/vol.tiff']);
    else
        vol = ref_volume_uneven(imio, opt, opt.nsamples, 17);  % caution: inline params
        % use one z slice for motion detection now;
        tiffwrite(vol, [opt.savepath '/vol.tiff']);
    end
    env.vol = vol;
    
    %% registration
    disp 'running registration';
    if opt.isregistration
%         if env.depth > 1
%             tic; rio = pipeline_registration_uneven(env, opt); toc;
%         else  % single core registration for single plane
            tic; rio = pipeline_registration_single_thread_gpu(env, opt); toc;
%         end
    else
        %
        temp = 1;
    end  % end if opt.isregistration
    
    %% f-b segmentation
    disp 'running background segmentation...'
    % caution: you may like to check the result of f-b segmentation
    [env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
    tiffwrite(im2uint8(env.volmask), [opt.savepath '/volmask.tiff']);
    tiffwrite(im2uint8(cat(2, rescalegd(env.vol, [0 0]), env.volmask)), [opt.savepath '/inspect_volmask.tiff']);
    %volmask = tiffread('volmask.tiff');
    %volmask = volmask > 100;
    
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
    
    %% save segmentation result
    %
    showspv = show_spv(env, opt, 'dot');
    svpath = checkpath([opt.savepath '/spvseg_dot/']);
    seqwrite(showspv, svpath);
    
    showspv = show_spv(env, opt, 'circle');
    svpath = checkpath([opt.savepath '/spvseg_circle/']);
    seqwrite(showspv, svpath);
    
    save([opt.savepath '/env.mat'], 'env', 'opt', '-v7.3');
    
    
    %% ---------------------------------------------------------------
    % extract signal
    disp 'extracting signals...';
    
    % % %     % this does not work because you must have a single imio in each
    % % %     % spmd thread and parellel pool does not support anonimous function
    % % %     % isregistration
    % % %     if opt.isregistration  % single plane
    % % %         readfunc = @(ti, zi) imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/' num2str(ti, '%.4d') '.tif']);
    % % %     else                   % no registration
    % % %         try  % it seems the file handler can expire
    % % %             close(imio);
    % % %         catch
    % % %         end
    % % %         imio = dcimgreader(opt.imgpath, opt.T, []);  %, opt.t0); % auto estimate t0
    % % %         readfunc = @(ti, zi) imio.readframe(ti, zi);
    % % %     end
    
    if opt.extract_signal_multicore  % extract with spmd
        
        %% extract with spmd
        nvolet = env.nvolet;
        height = env.height;
        width = env.width;
        depth = env.depth;
        
        mask = imdilate(env.mask, strel('disk', env.width/16+1));
        figure; imshow(mask)
        nspv = size(env.supervoxel, 1);
        [spvzind, spvregion] = supervoxel_index(env, opt);
        
        maskR = env.volmask;
        for dd=1:depth
            maskR(:,:,dd) = imdilate(env.volmask(:,:,dd), strel('disk', env.width/16+1));
            % figure(1000+dd); imshow(maskR(:,:,dd))
        end
        
        poolsize = 3;
        blocksize = nvolet; %ceil(nvolet / poolsize);
        
        % registered_io = env.registered_io;
        
        tic;
        %activities = zeros(nspv, env.nvol, 'single');
        %backgroundy = zeros(env.height, env.nvol, env.depth, 'single');
        
        poolobj = parpool(poolsize);
        
        if opt.isregistration
            
            spmd(poolsize)
                %registered_io = ImageIO([opt.savepath '/registered'], opt.T);
                %registered_io = dcimgreader('Z:\data\20151111\testdcimg.dcimg', opt.T, opt.t0);
                ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize nvolet]);
                
                actth = zeros(nspv, length(ti_in_this_thread));
                backth = zeros(height, length(ti_in_this_thread), depth);
                backthcf = zeros(length(ti_in_this_thread), depth);
                
                for ii = 1:length(ti_in_this_thread)
                    ti = ti_in_this_thread(ii);
                    [labindex ti]                   
                    
                    for zi = 1:depth
                        %         frame = dcimgread(registered_io, ti, zi);
                        if ~ismember(ti,opt.ind_stand_volet_total(:)) && ((zi<opt.ind_beg_step)||(zi>opt.ind_end_step))%%||(ismember(ti,234:252)&&ismember(zi,4:6))
                            backth(:,ii,zi) = nan;
                            backthcf(ii,zi) = nan;
                            for ppi = spvzind{zi}
                                actth(ppi, ii) = nan; %mean(frame(spvregion{ppi}));
                            end
                        else
                            frame = imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/f' num2str(ti, '%.4d') '.tiff']);
                            frame = single(frame);
                            frame = medfilt2(frame, [3 3]);  % caution: inline params; not sure whether should be here;
                        
                        % % %         if any(sum(~mask, 2) < 16)  % caution: inline params; not sure whether should be here;
                        % % %             backinty = zeros(size(mask,1), 1);
                        % % %         else
                        % % %             back = frame;
                        % % %             back(mask) = 0;
                        % % %             backinty = sum(back, 2) ./ (width  - sum(mask, 2));
                        % % %         end
                        
                            back = mean(frame(~mask));
                            if isnan(back)
                                back = 0;
                            end
                            backinty = back * ones(size(mask,1), 1);
                        
                            backth(:,ii,zi) = backinty;
                            %backgroundy(:, ti, zi) = backinty;
                        
                            backiz = mean(frame(~maskR(:,:,zi)));
                            if isnan(backiz)
                               backiz = 0;
                            end
                            backthcf(ii,zi) = backiz;
                        
                            % correct for each y
                            backinty = zeros(size(mask,1), 1);  % do not subtract the background
                            frame = frame - repmat(backinty, [1 width]);
                        
                            for ppi = spvzind{zi}
                                actth(ppi, ii) = mean(frame(spvregion{ppi}));
                            end
                        end
                    end % end zi
                    
                end % for ti
                if ~opt.isregistration
                    imio.close();
                end
            end  % spmd
            
        else  % ~opt.isregistration
            
            spmd(poolsize)
                
                % must have a reader in each spmd
                imio = dcimgreader(opt.imgpath, opt.T, []);  %, opt.t0); % auto estimate t0
                
                %registered_io = ImageIO([opt.savepath '/registered'], opt.T);
                %registered_io = dcimgreader('Z:\data\20151111\testdcimg.dcimg', opt.T, opt.t0);
                ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize nvolet]);
                
                actth = zeros(nspv, length(ti_in_this_thread));
                backth = zeros(height, length(ti_in_this_thread), depth);
                backthcf = zeros(length(ti_in_this_thread), depth);
                
                for ii = 1:length(ti_in_this_thread)
                    ti = ti_in_this_thread(ii);
                    [labindex ti]
                    
                    for zi = 1:depth
                        if ~ismember(ti,opt.ind_stand_vol_total(:)) && ((zi<opt.ind_beg_step)||(zi>opt.ind_end_step))
                            backth(:,ii,zi) = nan;
                            backthcf(ii,zi) = nan;
                            for ppi = spvzind{zi}
                                actth(ppi, ii) = mean(frame(spvregion{ppi}));
                            end
                        else
                        %         frame = dcimgread(registered_io, ti, zi);
                        %         frame = imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/' num2str(ti, '%.4d') '.tif']);
                        %         frame = readfunc(ti, zi);
                        frame = readframe(imio, opt, ti, zi);
                        frame = single(frame);
                        frame = medfilt2(frame, [3 3]);  % caution: inline params; not sure whether should be here;
                        
                        % % %         if any(sum(~mask, 2) < 16)  % caution: inline params; not sure whether should be here;
                        % % %             backinty = zeros(size(mask,1), 1);
                        % % %         else
                        % % %             back = frame;
                        % % %             back(mask) = 0;
                        % % %             backinty = sum(back, 2) ./ (width  - sum(mask, 2));
                        % % %         end
                        
                        back = mean(frame(~mask));
                        if isnan(back)
                            back = 0;
                        end
                        backinty = back * ones(size(mask,1), 1);
                        
                        backth(:,ii,zi) = backinty;
                        %backgroundy(:, ti, zi) = backinty;
                        
                        backiz = mean(frame(~maskR(:,:,zi)));
                        if isnan(backiz)
                            backiz = 0;
                        end
                        backthcf(ii,zi) = backiz;
                        
                        % correct for each y
                        backinty = zeros(size(mask,1), 1);  % do not subtract the background
                        frame = frame - repmat(backinty, [1 width]);
                        
                        for ppi = spvzind{zi}
                            actth(ppi, ii) = mean(frame(spvregion{ppi}));
                        end
                        end
                    end
                    
                end % for ti
                if ~opt.isregistration
                    imio.close();
                end
            end  % spmd
            
        end % if opt.isregistration
        duration_extract = toc;
        
        activities = cat(2, actth{:});
        backgroundy = cat(2, backth{:});
        backPSTH = cat(1,backthcf{:});
        save([opt.savepath '/act.mat'], 'activities', 'backgroundy','backPSTH', 'duration_extract', '-v7.3');
        
        delete(poolobj);
        
    else  % extract signals - single core
        
        %% extract signals - single core
        
        mask = imdilate(env.mask, strel('disk', env.width/16+1));
        
        tic;
        activities = zeros(size(env.supervoxel, 1), env.nvolet, 'single');
        backgroundy = zeros(env.height, env.nvolet, env.depth, 'single');
        for ti = 1:env.nvolet
            if mod(ti, 1000) == 1, disp(ti); end
            actt = zeros(size(env.supervoxel, 1), 1, 'single');
            for zi = 1:depth
                if ~ismember(ti,opt.ind_stand_vol_total(:)) && ((zi<opt.ind_beg_step)||(zi>opt.ind_end_step))
                    backgroundy (ti,zi) = nan;
                    for pi = spvzind{zi}
                        actt(pi) = nan;
                    end
                else
                %frame = imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/' num2str(ti, '%.4d') '.tif']);
                %         frame = readframe(imio, ti, zi);
                    frame = readfunc(ti, zi);
                    frame = single(frame);
                    frame = medfilt2(frame, [3 3]);  % caution: inline params; not sure whether should be here;
                
                    % % %         if any(sum(~mask, 2) < 16)  % caution: inline params; not sure whether should be here;
                    % % %             backinty = zeros(size(mask,1), 1);
                    % % %         else
                    % % %             back = frame;
                    % % %             back(mask) = 0;
                    % % %             backinty = sum(back, 2) ./ (width  - sum(mask, 2));
                    % % %         end
                
                    %         back = mean(frame(~mask));
                    %         if isnan(back)
                    %             back = 0;
                    %         end
                    back = 0;
                    backinty = back * ones(size(mask,1), 1);
                
                    %         back = frame;
                    %         back(mask) = 0;
                    %         backinty = sum(back, 2) ./ (width  - sum(mask, 2));
                
                    backgroundy(:, ti, zi) = backinty;
                
                    % correct for each y
                    frame = frame - repmat(backinty, [1 width]);
                
                   for pi = spvzind{zi}
                        actt(pi) = mean(frame(spvregion{pi}));
                   end
                end
            end
            activities(:,ti) = actt;
        end
        duration_extract = toc;
        save([opt.savepath '/act.mat'], 'activities', 'backgroundy', 'duration_extract', '-v7.3');
        
    end  % end if opt.extract_signal_multicore
    
%     batch code ------
end % batch for
% end batch code --




