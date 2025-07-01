    function [ROI_traces, ROIs] = test_ringfilter2_n
    % automatically chose ROIs
    % 201411_YangChen
    %20141201_cyin
    %20160831_cyin: for autoROI H2B-G6f nucleus expressing data; only differ in detect_by_ring_2d_n.m subfunction
%     cd('G:\Graduate\Academic\Graduate research\Du lab\Data_Imaging\20131219_SIN_Huc_G5\20131219_SIN_Huc_G5_exp_fish_5_38um')
% 	load('AVG_20131219_SIN_Huc_G5_exp_fish_5_38um_autoROI.mat')    
%%  parameters
    minrad = 3;
    maxrad = 8;
    thres = 0;

%%  
    fn_tif = uigetfile('*.tif'); % using tif stack file
%     stack
    stack = tiffread(fn_tif);
    % imread(fn_tif,index); tiffread.m same as imread.m for tiff stack frame by frame
    
    slice = mean(stack,3);
    slice = normim(slice);
    slice_back = slice;

    % smooth  smooth width should be smaller than cell diameter
    gs = fspecial('gaussian', minrad*2+1, 1);%
    slice = imfilter(slice, gs);

    slice = double(slice);
    show = slice / max(slice(:));

    radii = minrad:maxrad;
%     [points, fim, radmap] = detect_by_ring_2d_n(slice, radii, thres);
    [points, fim, radmap] = detect_by_ring_2d_n(slice, radii);    
    
    [y, x] = find(points);
    rad = sum(radmap, 3) / 2;
    rad = rad(points);
    
    rin = radmap(:,:,1);
    rin = rin(points);
    scores = fim(points);
    
    idx = remove_overlap([x y], rin, scores);
    
    rad = rad(idx);
    x = x(idx);
    y = y(idx);
    
    ind = sub2ind(size(slice), y, x);
    points = false(size(points));
    points(ind) = true;
    
    tshow = slice;
    tshow(points) = 1;
    
    figure; imshow(tshow);
    
    %{
    [y, x] = find(points);
    ind = find(points);
    rad = sum(radmap, 3) / 2;
    rad = rad(ind);
    %rad = trysize(diamap(ind))'/2-1;
    tshow = insertShape(tshow, 'circle', [x, y, rad], 'LineWidth', 1);
    %}

    %{
    conn = 4;
    cc = bwconncomp(points, conn);
    L = watershed_meyer(slice, conn, cc);

    show_ah = show;
    show_ah(L == 0) = 255;
    figure; imshow(show_ah);
    %}

    %%
    mask_region = roipoly; % manually define ROI area
    figure; imshow(mask_region);
    points(~mask_region) = false;
    
    [y, x] = find(points);
    rad = sum(radmap, 3) / 2;
    rad = rad(points);
    
    tshow = show;
    tshow(points) = true;

%     tshow = show;
%     tshow(points(:)) = 1;
    %{
    [y, x] = find(points);
    ind = find(points);
    rad = sum(radmap, 3) / 2;
    rad = rad(ind);
    %rad = trysize(diamap(ind))'/2-1;
    tshow = insertShape(tshow, 'circle', [x, y, rad], 'LineWidth', 1);
    %}
    
    f = figure('Name','autoROI','color', 'white', 'position', [400   200   size(slice)]);% 1 plot
    imshow(tshow);
    
%     text(x, y, num2str((1:length(x))'),'color','y');
    text(x, y, num2str((1:length(x))'));
	ellipse(rad,rad,zeros(size(rad)),x,y); 
    
    fn_tmp = fn_tif(1:end-15);
    fn = [fn_tmp '_autoROI_' num2str(length(x)) '_ROIs'];    
    saveas(f,fn,'tif')
    saveas(f,fn)
    
    for ii = 1:length(x)
        ROIs(ii).x = x(ii);
        ROIs(ii).y = y(ii);
        ROIs(ii).rad = rad(ii);
    end
% 	ind_discard = [];% manually discard false positive ROIs
%     ROIs(ind_discard) = [];
    
    %% time series of each ROI
    [xx, yy] = meshgrid(1:size(slice,2), 1:size(slice,1));
    ts = reshape(stack, [size(stack,1)*size(stack,2) size(stack,3)]);
    ROI_traces = ones(length(x), size(stack,3));
    for ii = 1:length(x)
        mask = (xx - x(ii)).^2 + (yy - y(ii)).^2 <= rad(ii)^2;
        ROI_traces(ii,:) = mean(ts(mask,:));
    end
% 	tt = imimposemin(show, points);% ·ÖË®ÁëËã·¨partition
%     L = watershed(tt, 8);
%     tt(L==0) = 1;
%     figure; imshow(tt)
            
    clear stack
    save(fn)
    end

            
