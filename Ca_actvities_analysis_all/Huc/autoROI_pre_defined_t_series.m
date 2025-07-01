    function ROI_traces = autoROI_pre_defined_t_series
    % automatically chose ROIs
    % 201411_YangChen
    %20141201_cyin
    %20150120_cyin: use previously determined ROIs to obtain time series
%%
    fn_mat = uigetfile('*.mat','previously autoROI mat file'); % using previously obtained autoROI mat file from test_ringfilter2.m
    load(fn_mat);

%%  
    fn_tif = uigetfile('*.tif'); % using tif stack file
%     stack
    stack = tiffread(fn_tif);
    % imread(fn_tif,index); tiffread.m same as imread.m for tiff stack frame by frame
    
    slice = mean(stack,3);

%     f = figure('Name','autoROI','color', 'white', 'position', [400   200   size(slice)]);% 1 plot
%     imshow(tshow);    
    
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

    fn_tmp = fn_tif(1:end-15);
    fn = [fn_tmp '_autoROI_' num2str(length(x)) '_ROIs'];  
    
    clear stack
%     save(fn)
    
    end

            
