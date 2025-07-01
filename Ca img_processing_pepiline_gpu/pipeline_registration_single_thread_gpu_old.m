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

dif_thres_z = zeros(depth, 1);
cor_thres_z = ones(depth, 1);
% 
warpopt = opt.warpopt;
%all_regtrans = cell(env.depth, env.nvolet); % save everytime

% todo: check whether reginfo exists, if it's not, create for all z
% poolsize = 12;
% poolobj = parpool(poolsize);

tic;
% spmd(poolsize)
%zi = 17;
imio = dcimgreader_uneven(opt.imgpath, opt);


for zi = 1:env.depth%%labindex:poolsize:env.depth  %zi = labindex + 12;
% [labindex zi]

zpath = checkpath([opt.savepath '/registered_new/z' num2str(zi, '%02d')]);
reginfopath = checkpath([opt.savepath '/registered_info/z' num2str(zi, '%02d')]);

cor_thres = cor_thres_z(zi);
dif_thres = dif_thres_z(zi);

% subsamplet = opt.subsamplet;
isreg = false(1, env.nvolet); 
ismove = true(1, env.nvolet);
dif_to_prev = zeros(1, env.nvolet);
cor_to_prev = ones(1, env.nvolet);

% continue
listdir = dir([reginfopath '/*.txt']);
t1 = 1;  %t1 = subsamplet(1);
if ~isempty(listdir)
    for iii = 1:length(listdir)
        filename = listdir(iii).name;
        fileid = str2double(filename(1:4));
        temp = dlmread([reginfopath  '/' filename]);
        isreg(fileid) = temp(1,1);
        ismove(fileid) = temp(1,2);
        cor_to_prev(fileid) = temp(2,1);
        dif_to_prev(fileid) = temp(2,2);
        t1 = max([t1 fileid]);  % subsamplet?
    end
end;

% ti = subsamplet(1);
refr_frame = env.vol(:,:,zi);
pre_refr = prepare_reference_image_for_registration_gpu(refr_frame, warpopt);
prev_frame = double(imio.readframe(opt, t1, zi));
[reg_frame, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, prev_frame, warpopt);
isreg(t1) = true;
prev_params = reg_params;

% save for ti=1
imwrite(uint16(reg_frame), [zpath '/' num2str(t1, '%04d') '.tif']);
dlmwrite([reginfopath  '/' num2str(t1, '%04d') '.txt'], ...
        [isreg(t1) isreg(t1); [cor_to_prev(t1) dif_to_prev(t1)]; squeeze(reg_params)]);
    
if (zi>=opt.ind_beg_step)&&(zi<=opt.ind_end_step)
    t2reg = t1:nvolet;
else
    t2reg = opt.ind_stand_volet_total(opt.ind_stand_volet_total>t1);
end

n_ti=length(t2reg);
for i_ti = 1:n_ti
    ti = t2reg(i_ti)  %subsamplet(2:end)  %2:4  %nvolet
    % compute similarity to last frame
    next_frame = double(imio.readframe(opt, ti, zi));
    dif_to_prev(ti) = mean((next_frame(:) - prev_frame(:)).^2);
    cor_to_prev(ti) = corr(next_frame(:), prev_frame(:));
    
    % decide whether to perform registration
    if cor_to_prev(ti) > cor_thres && dif_to_prev(ti) < dif_thres  % use transformation of last frame
        reg_frame = frame_reformat(next_frame, reg_params, warpopt);
    else  % register this frame
        [reg_frame, reg_params] = frame_registration_prepared_reference_gpu(pre_refr, next_frame, warpopt);
        isreg(ti) = true;
    end
    
    % ismove
    %if all(reg_params(:)==prev_params(:))  % see whether the frame is really moved, if not, update the thres
    if sum(abs(reg_params(:)-prev_params(:))) == 0 && ti > 50  % caution: inline params
        ismove(ti) = false;  % at least have one
    end
    
    % save registered image and registration information
    imwrite(uint16(reg_frame), [zpath '/' num2str(ti, '%04d') '.tif']);
    dlmwrite([reginfopath  '/' num2str(ti, '%04d') '.txt'], ...  % write every frame
        [[isreg(ti) ismove(ti)]; [cor_to_prev(ti) dif_to_prev(ti)]; [cor_thres dif_thres]; squeeze(reg_params)]);
    
    % prepare for next loop
    if ~ismove(ti)
        cor_thres = prctile(cor_to_prev(~ismove), prct);  % caution: inline params
        dif_thres = prctile(dif_to_prev(~ismove), 100-prct);  % one minutes
%         cor_thres = min(cor_to_prev(~ismove));  % caution: inline params
%         dif_thres = max(dif_to_prev(~ismove));  % one minutes
    end
    prev_frame = next_frame;
    prev_params = reg_params;
    %if ti > 100  % caution: inline params
    %    cor_thres = prctile(cor_to_prev(ti-100:ti-1), 40);
    %    dif_thres = prctile(dif_to_prev(ti-100:ti-1), 60);
    %end
end

end % for zi
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

