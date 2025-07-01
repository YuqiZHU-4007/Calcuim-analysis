function [duration_extract,activities,backgroundy,backPSTH ]=extract_signal_tiff(opt,env)
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
    nvol = env.nvol;
    height = env.height;
    width = env.width;
    depth = env.depth;
    
    mask = imdilate(env.mask, strel('disk', env.width/16+1));%%����mask
    h=figure('visible','off');imshow(mask);saveas(h,[opt.savepath '/mask2.tif']);
    nspv = size(env.supervoxel, 1);
    [spvzind, spvregion] = supervoxel_index(env, opt);
    
    maskR = env.volmask;
    for dd=1:depth
        maskR(:,:,dd) = imdilate(env.volmask(:,:,dd), strel('disk', env.width/16+1));
        % figure(1000+dd); imshow(maskR(:,:,dd))
    end
    
    poolsize = opt.poolsize;
    blocksize = ceil(nvol / poolsize);
    
    % registered_io = env.registered_io;
    
    tic;
    %activities = zeros(nspv, env.nvol, 'single');
    %backgroundy = zeros(env.height, env.nvol, env.depth, 'single');
    %delete(gcp('nocreate'))
    poolobj = parpool(poolsize);
    
    if opt.isregistration%%already registered above, load registered .tif files
        
        spmd(poolsize)
            %registered_io = ImageIO([opt.savepath '/registered'], opt.T);
            %registered_io = dcimgreader('Z:\data\20151111\testdcimg.dcimg', opt.T, opt.t0);
            ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize nvol]);
            
            actth = zeros(nspv, length(ti_in_this_thread));
            backth = zeros(height, length(ti_in_this_thread), depth);
            backthcf = zeros(length(ti_in_this_thread), depth);
            
            for ii = 1:length(ti_in_this_thread)
                ti = ti_in_this_thread(ii);
                [labindex ti]
                
                for zi = 1:depth
                    %         frame = dcimgread(registered_io, ti, zi);
                    frame = imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/' num2str(ti, '%.4d') '.tif']);
                    if length(find(frame~=0))==0
                        frame=nan(height,width);%, 'uint16'
                    else
                    frame = single(frame);
                    frame = medfilt2(frame, [3 3]);  % caution: inline params; not sure whether should be here;
                    end
                    % % %         if any(sum(~mask, 2) < 16)  % caution: inline params; not sure whether should be here;
                    % % %             backinty = zeros(size(mask,1), 1);
                    % % %         else
                    % % %             back = frame;
                    % % %             back(mask) = 0;
                    % % %             backinty = sum(back, 2) ./ (width  - sum(mask, 2));
                    % % %         end
                    
                    back = mean(frame(~mask));
                    if length(find(frame~=0))==0
                        backinty = back * ones(size(mask,1), 1);
                        backth(:,ii,zi) = backinty; 
                    else
                    if isnan(back)
                        back = 0;
                    end
                    backinty = back * ones(size(mask,1), 1);
                    
                    backth(:,ii,zi) = backinty;
                    end
                    %backgroundy(:, ti, zi) = backinty;
                    
                    backiz = mean(frame(~maskR(:,:,zi)));
                     if length(find(frame~=0))==0
                          backthcf(ii,zi) = backiz;
                     else
                    if isnan(backiz)
                        backiz = 0;
                    end
                    backthcf(ii,zi) = backiz;
                     end
                    % correct for each y
                    backinty = zeros(size(mask,1), 1);  % do not subtract the background, although calculated above
                    frame = frame - repmat(backinty, [1 width]);
                    
                    for ppi = spvzind{zi}
                        actth(ppi, ii) = mean(frame(spvregion{ppi}));
                    end
                end
                
            end % for ti
%             if ~opt.isregistration
%                 imio.close();
%             end
        end  % spmd
        
    else  % ~opt.isregistration %% do not register, load .dcimg file
        
        spmd(poolsize)
            
            % must have a reader in each spmd
            imio = dcimgreader_zyq_tiff_190827(opt.imgpath, opt.T, [],opt.process);  %, opt.t0); % auto estimate t0
            
            %registered_io = ImageIO([opt.savepath '/registered'], opt.T);
            %registered_io = dcimgreader('Z:\data\20151111\testdcimg.dcimg', opt.T, opt.t0);
            ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize nvol]);
            
            actth = zeros(nspv, length(ti_in_this_thread));
            backth = zeros(height, length(ti_in_this_thread), depth);
            backthcf = zeros(length(ti_in_this_thread), depth);
            
            for ii = 1:length(ti_in_this_thread)
                ti = ti_in_this_thread(ii);
                [labindex ti]
                
                for zi = 1:depth
                    %         frame = dcimgread(registered_io, ti, zi);
                    %         frame = imread([opt.savepath '/registered_new/z' num2str(zi, '%.2d') '/' num2str(ti, '%.4d') '.tif']);
                    %         frame = readfunc(ti, zi);
                    frame = readframe_zyq(imio, ti, zi);
                    if length(find(frame~=0))==0
                        actth(ppi, ii) = nan;
                    else
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
    %%δ�ģ�δ����֡��������
    
    mask = imdilate(env.mask, strel('disk', env.width/16+1));
    
    tic;
    activities = zeros(size(env.supervoxel, 1), env.nvol, 'single');
    backgroundy = zeros(env.height, env.nvol, env.depth, 'single');
    for ti = 1:env.nvol
        if mod(ti, 1000) == 1, disp(ti); end
        actt = zeros(size(env.supervoxel, 1), 1, 'single');
        for zi = 1:depth
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
        activities(:,ti) = actt;
    end
    duration_extract = toc; backPSTH=[];
    save([opt.savepath '/act.mat'], 'activities', 'backgroundy', 'duration_extract', '-v7.3');
    
end  % end if opt.extract_signal_multicore