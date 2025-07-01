function registered_io = pipeline_registration_single_thread_gpu(env, opt, prct)
% registered_io = pipeline_registration(env, opt)
%
% subsamplez and subsamplet will be ignored
% will not register those frames that they have similar appearance thresholded by
% estimating threshold from those frames being registered but do not change.
% 

%% test new registration
%dif_thres = 338.5078;
%cor_thres = 0.9832;

if nargin < 3
    prct = 5;
end

depth = env.depth;
nvolet = env.nvolet;

warpopt = opt.warpopt;

isreg = false(depth,nvolet); 
ismove = true(depth,nvolet);
dif_to_prev = zeros(depth, nvolet);
cor_to_prev = ones(depth, nvolet);
dif_thres = zeros(depth, nvolet);
cor_thres = ones(depth, nvolet);    
    
zpath = cell(env.depth,1);
reginfopath = cell(env.depth,1);

% poolsize = 12;
% poolobj = parpool(poolsize);

tic;
% spmd(poolsize)
%zi = 17;

imio = dcimgreader_uneven(opt.imgpath, opt);

refr_volume = gpuArray(single(env.vol(:,:,1:env.depth)));
pre_refr = prepare_reference_image_for_registration_gpu(refr_volume , warpopt);
prev_volume = single(nan(env.height,env.width,env.depth));

t1 = 1;  %t1 = subsamplet(1);
% todo: check whether reginfo exists, if it's not, create for all z
for zi = 1:env.depth%%labindex:poolsize:env.depth  %zi = labindex + 12;
    % [labindex zi]
    
    zpath{zi}=checkpath([opt.savepath '/registered_new/z' num2str(zi,'%02d') '/']);
    reginfopath{zi}=checkpath([opt.savepath '/registered_info/z' num2str(zi, '%02d') '/']);  
        
    % continue
    listdir = dir([reginfopath{zi} '/*.txt']);   
    if ~isempty(listdir)
            for iii = 1:length(listdir)
                filename = listdir(iii).name;
                fileid = str2double(filename(1:4));
                temp = dlmread([reginfopath{zi}  filename]);
                isreg(zi,fileid) = temp(1,1);
                ismove(zi,fileid) = temp(1,2);
                cor_to_prev(zi,fileid) = temp(2,1);
                dif_to_prev(zi,fileid) = temp(2,2);
                t1 = max([t1 fileid]);  % subsamplet?
            end            
    end
    
    prev_volume(:,:,zi) = single(imio.readframe(opt, t1, zi)); 

% ti = subsamplet(1);
end  %%zi
% the surface layer never targeted


prev_volume = gpuArray(single(prev_volume));
[reg_volume, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, prev_volume, warpopt); %%gpuArray in, gpuArray out
prev_params = reg_params;

isreg(:,t1) = true;

for zi = 1:env.depth
%     imwrite(uint16(gather(reg_volume(:,:,zi))),[opt.savepath '/registered_new/t' num2str(1,'%04d') '/z' num2str(zi,'%02d') '.tiff']);
%     t = Tiff([zpath{zi} '/f' num2str(t1,'%04d') '.tiff'],'w');
%     t.setTag('ImageLength',env.height);
%     t.setTag('ImageWidth', env.width);
%     t.setTag('Photometric', 1); %%MinIsBlack
%     t.setTag('BitsPerSample', 16);
%     t.setTag('SamplesPerPixel', 1);
%     t.setTag('TileWidth', 128);
%     t.setTag('TileLength', 128);
%     t.setTag('Compression', 1); %%None
%     t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%     t.setTag('Software', 'MATLAB');
%     t.write(uint16(gather(reg_volume(:,:,zi))));
%     t.close();
    imwrite(uint16(gather(reg_volume(:,:,zi))),[zpath{zi} '/f' num2str(t1,'%04d') '.tiff']);
    dlmwrite([reginfopath{zi}  '/' num2str(t1, '%04d') '.txt'], [isreg(zi,t1) ismove(zi,t1); [0 1]; [0 1]; squeeze(gather(reg_params(zi,:,:)))]);

end
% save for ti=1

for i_ti = find(opt.ind_stand_volet_total==t1)+1:length(opt.ind_stand_volet_total)
    next_volume = single(readvolume(imio, opt, opt.ind_stand_volet_total(i_ti)));
    next_volume = next_volume(:,:,1:env.depth);
    dif_to_prev_temp = gather(reshape(next_volume - prev_volume, env.height*env.width, env.depth));
    dif_to_prev_temp = mean(dif_to_prev_temp,1);
    dif_to_prev(:,opt.ind_stand_volet_total(i_ti)) = dif_to_prev_temp'.^2;
    for zi = 1:env.depth
        corrcoef_temp = corrcoef(next_volume(:,:,zi), prev_volume(:,:,zi));
        cor_to_prev(zi,opt.ind_stand_volet_total(i_ti)) = gather(corrcoef_temp(1,2));
    end
        
        % decide whether to perform registration
    if all(cor_to_prev(:,opt.ind_stand_volet_total(i_ti)) - cor_thres(:,opt.ind_stand_volet_total(i_ti))) && all(dif_to_prev(:,opt.ind_stand_volet_total(i_ti),:) < dif_thres(:,opt.ind_stand_volet_total(i_ti)))  % use transformation of last frame
       reg_volume = frame_reformat_gpu(next_volume, reg_params, warpopt);
    else  % register this frame
       [reg_volume, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, next_volume, warpopt);
       isreg(:,opt.ind_stand_volet_total(i_ti)) = true;
    end
        % ismove
        %if all(reg_params(:)==prev_params(:))  % see whether the frame is really moved, if not, update the thres
    if gather(sum(abs(reg_params(:)-prev_params(:)))) == 0 && i_ti > 50  % caution: inline params
       ismove(:,opt.ind_stand_volet_total(i_ti)) = false;  % at least have one
    end
     % save registered image and registration information     
    
   for zi = 1:env.depth
    imwrite(uint16(gather(reg_volume(:,:,zi))),[zpath{zi} '/f' num2str(opt.ind_stand_volet_total(i_ti),'%04d') '.tiff']);
%     t = Tiff([zpath{zi} '/f' num2str(opt.ind_stand_volet_total(i_ti),'%04d') '.tiff'],'w');
%     t.setTag('ImageLength',env.height);
%     t.setTag('ImageWidth', env.width);
%     t.setTag('Photometric', 1); %%MinIsBlack
%     t.setTag('BitsPerSample', 16);
%     t.setTag('SamplesPerPixel', 1);
%     t.setTag('TileWidth', 128);
%     t.setTag('TileLength', 128);
%     t.setTag('Compression', 1); %%None
%     t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%     t.setTag('Software', 'MATLAB');
%     t.write(uint16(gather(reg_volume(:,:,zi))));
%     t.close();
    dlmwrite([reginfopath{zi}  '/' num2str(opt.ind_stand_volet_total(i_ti), '%04d') '.txt'], ...  % write every frame
    [isreg(zi,opt.ind_stand_volet_total(i_ti)) ismove(zi,opt.ind_stand_volet_total(i_ti)); cor_to_prev(zi,opt.ind_stand_volet_total(i_ti)) dif_to_prev(zi,opt.ind_stand_volet_total(i_ti)); cor_thres(zi,opt.ind_stand_volet_total(i_ti)) dif_thres(zi,opt.ind_stand_volet_total(i_ti)); squeeze(gather(reg_params(zi,:,:)))]);    
    
   end
    
    for zi = 1:env.depth
    if ~ismove(opt.ind_stand_volet_total(i_ti))    
        indmove=find(~ismove(zi,:));
    % prepare for next loop   
        cor_thres(zi,opt.ind_stand_volet_total(min(i_ti+1,length(opt.ind_stand_volet_total)))) = prctile(cor_to_prev(zi,indmove), prct);  % caution: inline params
        dif_thres(zi,opt.ind_stand_volet_total(min(i_ti+1,length(opt.ind_stand_volet_total)))) = prctile(dif_to_prev(zi,indmove), 100-prct);  % one minutes
%         cor_thres = min(cor_to_prev(~ismove));  % caution: inline params
%         dif_thres = max(dif_to_prev(~ismove));  % one minutes
    end
    end
    
    prev_volume = next_volume;
    prev_params = reg_params;
    %if ti > 100  % caution: inline params
    %    cor_thres = prctile(cor_to_prev(ti-100:ti-1), 40);
    %    dif_thres = prctile(dif_to_prev(ti-100:ti-1), 60);
    %end
end

for iTrial = 1:opt.nTrial
        
    refr_volume = gpuArray(single(env.vol(:,:,opt.ind_beg_step:opt.ind_end_step)));
    pre_refr = prepare_reference_image_for_registration_gpu(refr_volume , warpopt);
    
    ind_last_volet = (iTrial-1)*opt.nvolets_per_trial+opt.ind_og_vol;
    prev_volume = readvolume(imio, opt, ind_last_volet);    
    prev_volume = gpuArray(single(prev_volume(:,:,opt.ind_beg_step:opt.ind_end_step)));
    
    [~, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, prev_volume, warpopt); %%gpuArray in, gpuArray out
    prev_params = reg_params;

    for i_ti = 1:opt.nsteps_og_vol/(opt.ind_end_step-opt.ind_beg_step)/2
        next_volume = readvolume(imio, opt, ind_last_volet+i_ti);
        next_volume = gpuArray(single(next_volume(:,:,opt.ind_beg_step:opt.ind_end_step)));
        dif_to_prev_temp = gather(reshape(next_volume - prev_volume, env.height*env.width, opt.ind_end_step-opt.ind_beg_step+1));
        dif_to_prev_temp = mean(dif_to_prev_temp,1);
        dif_to_prev(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti) = dif_to_prev_temp'.^2;
        for zi = 1:opt.ind_end_step-opt.ind_beg_step+1
            corrcoef_temp = corrcoef(next_volume(:,:,zi), prev_volume(:,:,zi));
            cor_to_prev(zi,ind_last_volet+i_ti) = gather(corrcoef_temp(1,2));
        end
                
        % decide whether to perform registration
        if all(cor_to_prev(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti) - cor_thres(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti)) && all(dif_to_prev(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti,:) < dif_thres(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti))  
            % use transformation of last frame
            reg_volume = frame_reformat_gpu(next_volume, reg_params, warpopt);
        else
            % register this frame
            [reg_volume, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, next_volume, warpopt);
            isreg(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti) = true;
        end
        % ismove
        %if all(reg_params(:)==prev_params(:))  % see whether the frame is really moved, if not, update the thres
        if gather(sum(abs(reg_params(:)-prev_params(:)))) == 0 && i_ti > 50  % caution: inline params
           ismove(opt.ind_beg_step:opt.ind_end_step,ind_last_volet+i_ti) = false;  % at least have one
        end
        % save registered image and registration information
        
        for zi = 1:opt.ind_end_step-opt.ind_beg_step+1
            imwrite(uint16(gather(reg_volume(:,:,zi))),[zpath{opt.ind_beg_step+zi-1} '/f' num2str(ind_last_volet+i_ti,'%04d') '.tiff']);
%     t = Tiff([zpath{zi} '/f' num2str(opt.ind_stand_volet_total(i_ti),'%04d') '.tiff'],'w');
%     t.setTag('ImageLength',env.height);
%     t.setTag('ImageWidth', env.width);
%     t.setTag('Photometric', 1); %%MinIsBlack
%     t.setTag('BitsPerSample', 16);
%     t.setTag('SamplesPerPixel', 1);
%     t.setTag('TileWidth', 128);
%     t.setTag('TileLength', 128);
%     t.setTag('Compression', 1); %%None
%     t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%     t.setTag('Software', 'MATLAB');
%     t.write(uint16(gather(reg_volume(:,:,zi))));
%     t.close();
            dlmwrite([reginfopath{opt.ind_beg_step+zi-1}  '/' num2str(ind_last_volet+i_ti, '%04d') '.txt'], ...  % write every frame
    [isreg(zi,ind_last_volet+i_ti) ismove(zi,ind_last_volet+i_ti); cor_to_prev(zi,ind_last_volet+i_ti) dif_to_prev(zi,ind_last_volet+i_ti); cor_thres(zi,ind_last_volet+i_ti) dif_thres(zi,ind_last_volet+i_ti); squeeze(gather(reg_params(zi,:,:)))]);    
        end
        
       for zi = opt.ind_beg_step:opt.ind_end_step
           if ~ismove(ind_last_volet+i_ti)  
               indmove=find(~ismove(zi,:));
               % prepare for next loop   
               cor_thres(zi,ind_last_volet+i_ti+1)= prctile(cor_to_prev(zi,indmove), prct);  % caution: inline params
               dif_thres(zi,ind_last_volet+i_ti+1) = prctile(dif_to_prev(zi,indmove), 100-prct);  % one minutes
               %         cor_thres = min(cor_to_prev(~ismove));  % caution: inline params
               %         dif_thres = max(dif_to_prev(~ismove));  % one minutes
           end
       end
    
    prev_volume = next_volume;
    prev_params = reg_params;
    %if ti > 100  % caution: inline params
    %    cor_thres = prctile(cor_to_prev(ti-100:ti-1), 40);
    %    dif_thres = prctile(dif_to_prev(ti-100:ti-1), 60);
    end
end

imio.close();
% end % spmd
%%---------------------------------------------


%{
%% old registration 
zi = 1; ti = 1;
parfor ti = 1:24
    ti
    frame = readframe(imio, ti, zi);
end

regscores = zeros(env.depth, env.nvolet);
regtrans  = cell(env.depth, env.nvolet);
for zi = opt.subsamplez
    zi
    tic;
    zpath = checkpath([opt.savepath 'registered/z' num2str(zi, '%02d')]);
    parfor ti = 1367:env.nvolet  %     for ti = 1:env.nvolet
        ti
        %if exist([zpath '/' num2str(ti, '%04d') '.tif'], 'file') > 0
        %    continue;
        %end
        frame = readframe(env.imio, ti, zi);
        %frame = imior.readframe(ti, zi);
        [mvreg, trans] = frame_registration(frame, env.vol(:,:,zi), warpopt);
        diffim = mvreg - env.vol(:,:,zi);
        regscores(zi, ti) = mean(diffim(:));
        regtrans{zi, ti} = trans;
        imwrite(mvreg, [zpath '/' num2str(ti, '%04d') '.tif']);
    end
    toc;
end  % for zi
%}

% delete(poolobj);
registered_io = [];%ImageIO([opt.savepath '/registered_new'], opt.T); 

