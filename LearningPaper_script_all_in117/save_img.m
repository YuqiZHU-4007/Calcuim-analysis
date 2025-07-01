function save_img(img,filename,savepath)
pictureDir = savepath;
% 定义一个结构体，用于放生成动画的帧
movie_i = struct;
movie_i.cdata = [];
movie_i.colormap = [];

%% 逐个读取图片信息，并保存动画的帧
% 设置一个figure，用于显示动画：根据自己的图片尺寸进行调节
fig = figure('position',[100,100,size(img,2),size(img,1)]);
len = size(img,4);
% 逐个读取图片
for iP=1:len
    picture = img(:,:,:,iP); % 导入图片
    picture_fy3_resize = imresize(picture,1); % 对图片进行重采样
    movie_i(iP) = im2frame(picture_fy3_resize); % 图片保存为动画的帧
end
% 在fig中播放动画，播放1遍，速度1帧/秒
movie(fig,movie_i ,1,5)

%% 生成视频，并保存
 writerObj =VideoWriter(fullfile(savepath,filename)); % 生成一个avi动画
 writerObj.FrameRate=5; % 设置avi动画的参数，设置帧速率
 open(writerObj); % 打开avi动画
 writeVideo(writerObj,movie_i); % 将保存的动画写入到视频文件中
 close(writerObj); % 关闭动画
 close(fig);
%  
%  a=img(:,:,:,30);
%  figure,imshow(a,[min(a(:)) max(a(:))/10])