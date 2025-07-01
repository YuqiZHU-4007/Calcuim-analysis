  clc;clear;close all;

    fn_tif = uigetfile('*.tif');
    ori_stack = tiffread(fn_tif,'uint16');
%     stack = im2double(ori_stack);
    [y,x,z] = size(ori_stack); % 纵轴为y,行
    imagesc(ori_stack(:,:,1));
    roi_size =0.5; % micron
    pixel_size = 0.1036; % micron
    bin = round(roi_size/pixel_size);

    iy = ceil(y/bin)-1;
    ix = ceil(x/bin)-1;
    seg_stack = zeros(iy,ix,z);
    dFvsFo = zeros(iy,ix,z);
    smoothwin = 100; 
    num_of_std = 3;
    frame_time = 0.061597;% in s
    tmp_F = zeros(z,iy*ix);
    dFvsFo_all = zeros(z,iy*ix);
    for i = 1:iy 
        for j = 1:ix
            for t = 1:z    
                 tmp = ori_stack((i-1)*bin+1:i*bin,(j-1)*bin+1:j*bin,t);
                 tmp_F(t,(i-1)*ix+j) = sum(sum(tmp))/bin^2;     %算出每个superpixel随时间变化的F值
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
       
% F_t = cell (iy,ix); % 每个superpixel每一帧的值
% for i = 1:iy
%     for j = 1:ix
%         F_t{i,j} = tmp_F(:,(i-1)*ix + j);
%     end
% end       

            
    %%
    all_F = sum(sum(ori_stack,1),2)/(x*y);
    tmp = [];
    h = figure();
    for t = 1:z
        subplot(2,1,1);
        tmp = [tmp,all_F(t)];
        time = frame_time:frame_time:t*frame_time;
        plot(time,tmp);
        xlim([0 z*frame_time]);
        ylim([min(all_F),max(all_F)]),
        xlabel('time (s)');
        ylabel('F');
         hold on 
%      plot([41,41],[200,240],'r'); %looming
%      plot([137,137],[200,240],'r');
%      plot([233,233],[200,240],'r');
%      plot([329,329],[200,240],'r');
%      plot([425,425],[200,240],'r');

%        plot([59,59],[153,160],'r'); % food extract
%        plot([119,119],[153,160],'r');
%        plot([179,179],[153,160],'r');
%        plot([239,239],[153,160],'r');
%        plot([299,299],[153,160],'r');
        hold off
        box off

        subplot(2,1,2);
        imshow(seg_stack(:,:,t),[0 0.5]);
%         imshow(seg_stack(:,:,t),[0 1]); % 匹配norm_dFvsFo
        colormap(hot);
        colorbar;
        F(t) = getframe(h); % online show
        %saveas(h,[num2str(t) 'bin=4'],'tif');
        
%         imshow(seg_stack(:,:,t),[0 0.5]);
%         saveas(gcf,num2str(t),'tif');
%         close all
%         disp(num2str(t));
    end 
    obj = VideoWriter('whole-brain_bin=5.avi','Uncompressed AVI');
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

% peak_onset = cell(1,1);  %根据全局的F值找出DA release event onset
% sliding_bin = 20; % sliding window, in frame
% for i =1
%     peak_onset{i} = [];
%     tmp = all_F;
%     for j = 1:z-2*sliding_bin-1
%         baseline_mean = mean(tmp(j:j+sliding_bin-1)); % j+sliding_bin is the checked frame
%         baseline_std = std(tmp(j:j+sliding_bin-1));
%         if tmp(j+sliding_bin) > baseline_mean + 2 * baseline_std && tmp(j+sliding_bin+1) > baseline_mean + 2 * baseline_std && tmp(j+sliding_bin) > 4
%            peak_onset{i} = [peak_onset{i},j+sliding_bin];
%         end
%     end
%      tmptmp = peak_onset{i} - [peak_onset{i}(1),peak_onset{i}(1:end-1)];
%      onset_id = find(tmptmp~=1);
%      peak_onset{i} = peak_onset{i}(onset_id);
% end
% 
%  %根据全局的F值找出DA release event结束的帧数
% baseline_mean = []; 
% peak_duration = [];
% peak_off = [];
% for i = 1:length(peak_onset{1})
%     baseline_mean(i) = mean(all_F(1,1,peak_onset{1}(i)-sliding_bin:peak_onset{1}(i)-1));
%     baseline_std(i) = std(all_F(1,1,peak_onset{1}(i)-sliding_bin:peak_onset{1}(i)-1));
%     for k = 1:2*sliding_bin
%         if all_F(1,1,peak_onset{1}(i)+k) > baseline_mean(i) + 2 * baseline_std(i) && all_F(1,1,peak_onset{1}(i)+k+1) < baseline_mean(i) + 2 * baseline_std(i)  
%         peak_duration(i) = k+1;
%         peak_off(i) = peak_onset{1}(i)+k+1;
%         end
%     end
% end
% 
% peak_on_off = cell(1,length(peak_onset{1}));
% for i = 1:length(peak_onset{1})
%     peak_on_off{i} = [peak_onset{1}(i),peak_off(i)];
% end
%  

% peakidentification_bin = 20;
% tmp = all_F;
% for t = 1:z - peakidentification_bin - 1
%     baseline_t = mean(tmp(t:t+peakidentification_bin-1));
%     baseline_std = std(tmp(t:t+peakidentification_bin-1));
%     if tmp(t+peakidentification_bin) > baseline_t + 3 * baseline_std && tmp (t+peakidentification_bin+1) > baseline_t + 3 * baseline_std && tmp(t+peakidentification_bin) > 4
%        peak_onset =[peak_onset, j + peakidentification_bin];         
%     end
%     tmptmp = peak_onset - [peak_onset(1),peak_onset(1:end-1)];
%     onset_id = find(tmptmp~=1);
%     peak_onset = peak_onset(onset_id);
% end
%     tmptmp = peak_onset - [peak_onset(1),peak_onset(1:end-1)];
%     onset_id = find(tmptmp~=1);
%     peak_onset = peak_onset(onset_id);             
        
     
% 识别反应onset，按peak提取每个release event的spatial pattern
%     peakidentification_bin = 20;
%     peak_onset = cell(iy,ix);
%     for i = 1:iy
%         for j = 1:ix
%             peak_onset{i,j} = [];
%             tmp = smooth(detrend(F_t{i,j}),3);
%             % kk=1;
%             for t = 1:z - peakidentification_bin -1
%             baseline_t = mean(tmp(t:t+peakidentification_bin-1));
%             baseline_std = std(tmp(t:t+peakidentification_bin-1));
%             if tmp(t+peakidentification_bin) > baseline_t + 3 * baseline_std && tmp (t+peakidentification_bin+1) > baseline_t + 3 * baseline_std && tmp(t+peakidentification_bin) > 4
%                 % peak_onset{iy,ix} (kk)=[t+peakidentification_bin];
%                 peak_onset{iy,ix} =[peak_onset{iy,ix}, j + peakidentification_bin];
%                 % kk=kk+1;
%             end
%             end
% %              tmptmp = peak_onset{iy,ix} - [peak_onset{iy,ix}(1),peak_onset{iy,ix}(1:end-1)];
% %              onset_id = find(tmptmp~=1);
% %              peak_onset{i,j} = peak_onset{iy,ix}(onset_id)
%          end
%     end
     
    
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


        
% figure,
% dFvsFo_all_corr =dFvsFo_all(:,:);
% [R,P] = corrcoef(dFvsFo_all_corr);
% P = 1 - P;
% clims = [-0.5 0.5];
% f = imagesc(R,clims);

% 每个superpixel每一帧的dFvsFo
dFvsFo_t = cell (iy, ix); 
for i = 1:iy 
    for j = 1:ix
          dFvsFo_t{i, j} = dFvsFo_all (:,(i-1)*ix +j);
    end
end

num_currentinjection = 5;
dFvsFo_currentinjection = struct(1,num_currentinjection);
for i = 1:num_currentinjection
end




% 画出每个peak对应有反应的superpixel pattern
peak_dFvsFo_t = cell(1,length(peak_onset{1}));
mean_dFvsFo_t_peak = cell(1,length(peak_onset{1}));
for k = 1:length(peak_onset{1})
    peak_dFvsFo_t{k} = dFvsFo_all(peak_onset{1}(k):peak_off(k),:);
    mean_dFvsFo_t_peak{k} = zeros(iy,ix);
    for i = 1:iy    
        for j = 1:ix
            dFvsFo_t_peaks = peak_dFvsFo_t{k}(:,(i-1)*ix+j);
%             mean_dFvsFo_t_peak{k}(i,j) = sum(dFvsFo_t_peaks)/(peak_off(k)-peak_onset{1}(k)+1);
            mean_dFvsFo_t_peak{k}(i,j) = mean(dFvsFo_t_peaks);
        end
    end
    figure, imagesc(mean_dFvsFo_t_peak{k},[-0.01 0.1]);
    %set gca "position" [x,y,10,20]
end

%自定义画出某些帧的平均
peak_dFvsFo_t1 = [];
peak_dFvsFo_t1 = dFvsFo_all(1000:2000,:);
dFvsFo_t_peak = cell (iy, ix); % 每个superpixel每个peak中每一帧的dFvsFo
mean_dFvsFo_t_peak = zeros (iy,ix);
for i = 1:iy 
    for j = 1:ix
               dFvsFo_t_peak{i, j} = peak_dFvsFo_t1 (:,(i-1)*ix +j);
               mean_dFvsFo_t_peak(i,j) = sum(dFvsFo_t_peak{i, j})/34;
    end
end
figure, imagesc(mean_dFvsFo_t_peak,[-0.01 0.05]);







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
mean_dFvsFo_eachROI = zeros(iy,ix);
for i = 1:iy, 
    for j = 1:ix
        mean_F_eachROI(i,j) = sum(tmp_F(:,(i-1)*ix+j))/z;
    end
end

% 每个superpixel的dFvsFo均值
mean_dFvsFo_eachROI = zeros(iy,ix);
for i = 1:iy, 
    for j = 1:ix
        mean_dFvsFo_eachROI(i,j) = sum(dFvsFo_all(:,(i-1)*ix+j))/z;
    end
end

figure,  scatter(reshape(mean_F_eachROI,[],1),reshape(mean_dFvsFo_eachROI,[],1));
for i = 1:iy
    for j = 1:ix
        scatter(mean_F_eachROI(i,j),mean_dFvsFo_eachROI(i,j));
        hold on
    end
end



figure, surf(1:ix,1:iy,mean_dFvsFo_eachROI);
zlim([-0.001,0.001]);

divide_dFvsFo_F = zeros(iy,ix);
for i = 1:iy,
    for j = 1:ix
        divide_dFvsFo_F(i,j) = mean_dFvsFo_eachROI(i,j)/mean_dFvsFo_eachROI(i,j);
    end
end

h=figure, imagesc(mean_dFvsFo_eachROI,[0 0.005]);
saveas(h,'mean_dFvsFo_eachROI','pgm');

k = figure, imagesc(mean_F_eachROI,[80 700]);
saveas(k,'mean_F_eachROI','pgm');

figure, imagesc(mean_dFvsFo_eachROI./mean_F_eachROI,[0 0.0001]);



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
A = imrotate(dFvsFo_sig_deletelowmeanvalue, -12);
figure, imagesc(A,[0 1]);
xlim([0 80]);
ylim([5 140]);







    