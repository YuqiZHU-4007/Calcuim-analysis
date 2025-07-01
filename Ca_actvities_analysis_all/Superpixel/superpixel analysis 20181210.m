  clc;clear;close all;

    fn_tif = uigetfile('*.tif');
    ori_stack = tiffread(fn_tif,'uint16');
%     stack = im2double(ori_stack);
    [y,x,z] = size(ori_stack); % 纵轴为y,行
    imagesc(ori_stack(:,:,1));
    roi_size = 0.397; % micron
    pixel_size = 0.397; % micron
    bin = round(roi_size/pixel_size);

    iy = ceil(y/bin)-1;
    ix = ceil(x/bin)-1;
    seg_stack = zeros(iy,ix,z);
    dFvsFo = zeros(iy,ix,z);
    smoothwin = 100;
    num_of_std = 3;
    frame_time = 1;
    tmp_F = zeros(z,iy*ix);
    dFvsFo_all = zeros(z,iy*ix);
    for i = 1:iy 
        for j = 1:ix
            for t = 1:z    
                 tmp = ori_stack((i-1)*bin+1:i*bin,(j-1)*bin+1:j*bin,t);
                 tmp_F(t,(i-1)*ix+j) = sum(sum(tmp))/bin^2;     
%                  seg_stack(i,j,t) = sum(sum(tmp))/bin^2;
            end
            data_this_ROI = smooth(tmp_F(:,(i-1)*ix+j),3);
            smoothfactor = round(smoothwin / frame_time);
            baseline_ind_iteration = 1:length(data_this_ROI);  % Logical vector holding the indices of baseline for each iteration
            baseline_discard = 0;
            count = 0;     
            while ~isempty(baseline_discard) && count <= 100
                  baseline_smoothed = smooth(data_this_ROI(baseline_ind_iteration),smoothfactor);
                  baseline_smoothed_std = std(baseline_smoothed - data_this_ROI(baseline_ind_iteration));
                  baseline_threshold = baseline_smoothed + baseline_smoothed_std * num_of_std; % mean + n-fold * sd
                  baseline_discard = find(data_this_ROI(baseline_ind_iteration) > baseline_threshold); % 这些帧不作为baseline？
                  baseline_ind_iteration(baseline_discard) = []; % 去除了部分帧
                  peak_ind = setdiff(1:length(data_this_ROI),baseline_ind_iteration); % 疑似反应对应的帧
                  count = count + 1;
            end
            baseline_index = unique([1 baseline_ind_iteration length(data_this_ROI)]); % baseline对应的帧
            data_smooth_baseline = interp1(baseline_index,smooth(baseline_index,data_this_ROI(baseline_index),smoothfactor),1:length(data_this_ROI));
            data_baseline_removed = (data_this_ROI - data_smooth_baseline')./data_smooth_baseline'; % 重要
            dFvsFo = data_baseline_removed; % 重要
            dFvsFo_all(:,(i-1)*ix+j) = dFvsFo;
            for tt = 1:z
                seg_stack(i,j,tt) = dFvsFo(tt);
            end
%             norm_dFvsFo = dFvsFo/max(dFvsFo); % 匹配norm_dFvsFo
%             dFvsFo_all(:,(i-1)*ix+j) = norm_dFvsFo; % 匹配norm_dFvsFo
%             for tt = 1:z % 匹配norm_dFvsFo
%                 seg_stack(i,j,tt) = norm_dFvsFo(tt); % 匹配norm_dFvsFo
%             end % 匹配norm_dFvsFo
            disp([num2str(i),' ',num2str(j)]);
        end
    end
    
    %%
    all_F = sum(sum(ori_stack,1),2)/(x*y);
    tmp = [];
    h = figure();
    for t = 1:z
        subplot(2,1,1);
        tmp = [tmp,all_F(t)];
        time = 1:1:t;
        plot(time,tmp);
        xlim([0 z]);
        ylim([min(all_F),max(all_F)]),
        xlabel('time ( frame )');
        ylabel('F');
        box off

        subplot(2,1,2);
        imshow(seg_stack(:,:,t),[0 0.5]);
%         imshow(seg_stack(:,:,t),[0 1]); % 匹配norm_dFvsFo
        colormap(hot);
        colorbar;
        F(t) = getframe(h); % online show
        
%         imshow(seg_stack(:,:,t),[0 0.5]);
%         saveas(gcf,num2str(t),'tif');
%         close all
%         disp(num2str(t));
    end 
    obj = VideoWriter('whole-brain_bin=20.avi','Uncompressed AVI');
    open(obj)
    writeVideo(obj,F);
    close(obj);
%     h = figure();
% for t = 1:z
%     imagesc(seg_stack(:,:,t),[0 0.5]);
%     colormap(hot);
%     colorbar;
%     F(t) = getframe(h);
% end
   
%%
%     h = figure();
%     v = VideoWriter('myIndexed.avi','Indexed AVI');
% for k = 1:50
%     imagesc(seg_stack(:,:,k));
%     colormap(jet); 
%     saveas(gcf,num2str(k),'tif');
%     close all
%     A = imread([num2str(k),'.tif']);
%     [X,map] = rgb2ind(A,64);
%     v.Colormap = map;
%     open(v)
%     writeVideo(v,X)
%     close(v)
% end
% 
% obj = VideoWriter('1.avi','Uncompressed AVI');
% open(obj)
% writeVideo(obj,rgb);
% close(obj);
%%
    
%     F = zeros(z,ix * iy);
%     max_dFvsF0_all = zeros(iy,ix);
%     for i = 1:iy
%         for j = 1:ix
%             for t = 1:z
%                 F(t,(i - 1) * j + j) = seg_stack(i,j,t);
%             end
%             tmp = F(:,(i - 1) * j + j);
%             baseline = sort(tmp,'ascend');
%             baseline = mean(baseline(1:10));
%             max_dFvsF0 = (max(tmp) - baseline)/baseline;
%             max_dFvsF0_all(i,j) = max_dFvsF0;
%         end
%     end
    
%     [frame_num,num_ROI] = size(F);
%     for i = 1:num_ROI
%         tmp = F(:,i);
%         
%     end

%     imagesc(max_dFvsF0_all,[0,2]);
%     fn_tif_marker = uigetfile('*.tif');
%     ori_stack = tiffread(fn_tif_marker);
%     figure
%     a = mean(ori_stack,3);
%     a(a<180) = 0;
%     imagesc(a,[0,400]);
%     b = max_dFvsF0_all*0.5+a*0.5;
%     imshow(b);
    
%     clim = [0,10];
%     imagesc(a',clim);
%     imagesc(dFvsFo');
%     for i = 1:47
%         figure
%     a = seg_stack(i,2,:);
%     plot(smooth(a,20));
%     end

dFvsFo_t = cell (iy, ix); % 每个superpixel每一帧的值
for i = 1:iy 
    for j = 1:ix
               dFvsFo_t{i, j} = dFvsFo_all (:,(i-1)*ix +j);
    end
end
        
figure,
dFvsFo_all_corr =dFvsFo_all(:,:);
[R,P] = corrcoef(dFvsFo_all_corr);
P = 1 - P;
clims = [-0.5 0.5];
f = imagesc(R,clims);

reshape_dFvsFo_all = dFvsFo_all(:); % 所有像素点所有时间点的dFvsFo捋成一列数组
% deleteNaN_reshape_dFvsFo_all = rmmissing(reshape_dFvsFo_all);

deletenan_reshape_dFvsFo_all = reshape_dFvsFo_all(isfinite(reshape_dFvsFo_all)); % 去除NaN的点

%  norm_deletenan_reshape_dFvsFo_all = normfit(deletenan_reshape_dFvsFo_all);
%  [miu,sigma] = normfit(deletenan_reshape_dFvsFo_all);
%  threshold_sig = miu + sigma;

sort_deletenan_reshape_dFvsFo_all = sort(deletenan_reshape_dFvsFo_all); % 所有superpixel所有dFvsFo值从小到大排列
A = sort_deletenan_reshape_dFvsFo_all;
threshold_sig = A(round(0.8 * length(A))); % 找到前20%的dFvsFo

% 每个superpixel的荧光亮度均值
mean_F_eachROI = zeros(iy,ix);
for i = 1:iy, 
    for j = 1:ix
        mean_F_eachROI(i,j) = sum(tmp_F(:,(i-1)*ix+j))/z;
    end
end
% 根据ROI的基底荧光亮度值筛选，过低的刨除
reshape_mean_F_eachROI = mean_F_eachROI(:);
sort_reshape_mean_F_eachROI = sort(reshape_mean_F_eachROI);
% threshold_sort_reshape_mean_F_eachROI = sort_reshape_mean_F_eachROI(round(0.8*length(sort_reshape_mean_F_eachROI)));
threshold_minimalmeanvalue = 80; %原始数据的mean值至少为100才被纳入

% 每个superpixel的dFvsFo大于设定阈值的帧数
dFvsFo_sig = zeros(iy,ix);
for i = 1:iy
    for j = 1:ix
        tmp = dFvsFo_t{i,j};     
        num = length(find(tmp>threshold_sig));
        dFvsFo_sig(i,j) = num;
    end
end

% 去除基底荧光亮度过低的superpixel
dFvsFo_sig_deletelowmeanvalue = zeros (iy,ix);
for i = 1:iy
    for j = 1:ix
        if mean_F_eachROI(i,j) < threshold_minimalmeanvalue;
            dFvsFo_sig_deletelowmeanvalue(i,j) = 0;
        else dFvsFo_sig_deletelowmeanvalue (i,j) = dFvsFo_sig(i,j);
        end
    end
end
 
dFvsFo_sig = dFvsFo_sig/max(max(dFvsFo_sig));
imagesc(dFvsFo_sig,[0 1]);

figure, imagesc(dFvsFo_sig,[0 1]);


dFvsFo_sig_deletelowmeanvalue = dFvsFo_sig_deletelowmeanvalue/max(max(dFvsFo_sig_deletelowmeanvalue));
figure, imagesc(dFvsFo_sig_deletelowmeanvalue,[0 1]);
A = imrotate(dFvsFo_sig_deletelowmeanvalue, 15);
figure, imagesc(A,[0 1]);
xlim([5 100]);
ylim([10 130]);






A1 = dFvsFo_t{20,11};
A2 = dFvsFo_t{20,14};
A3 = dFvsFo_t{21,19};
A4 = dFvsFo_t{50,18};

A = [A1 A2 A3 A4];
xlswrite('dFvsFofoursuperpixel.xlsx',A);

% all_F_t = zeros(z,y*x);
% for i = 1:z
%     for j = 1:y*x
%         if j ~= idivide (y*x, x, 'floor')
%             all_F_t(i,j) = ori_stack(floor(j/504)+1,j-floor(j/504)*504,i);
%         else all_F_t(i,j) = ori_stack(floor(j/504),j-(floor(j/504)-1)*504,i);
%         end
%     end
% end

%每个像素点（pixel)随时间变化的亮度值
all_F_t = zeros(y,x,t);
for i = 1:y
    for j = 1:x
        for t = 1:z
        all_F_t(i,j,t) = ori_stack(i,j,t);
        end
    end
end

F_everysuperpixel2011 = zeros (z,bin^2);
for i = 1:z
    for j = 1:bin^2
        F_everysuperpixel2011(i,j) = all_F_t(201+(ceil(j/bin)-1),110+j-(ceil(j/bin)-1)*bin,i);
    end
end

xlswrite('F_everysuperpixel2011.xlsx',F_everysuperpixel2011);
        
F_everysuperpixel2021 = zeros (z,bin^2);
for i = 1:z
    for j = 1:bin^2
        F_everysuperpixel2021(i,j) = all_F_t(201+(ceil(j/bin)-1),210+j-(ceil(j/bin)-1)*bin,i);
    end
end

xlswrite('F_everysuperpixel2021.xlsx',F_everysuperpixel2021);



    