function save_img(img,filename,savepath)
pictureDir = savepath;
% ����һ���ṹ�壬���ڷ����ɶ�����֡
movie_i = struct;
movie_i.cdata = [];
movie_i.colormap = [];

%% �����ȡͼƬ��Ϣ�������涯����֡
% ����һ��figure��������ʾ�����������Լ���ͼƬ�ߴ���е���
fig = figure('position',[100,100,size(img,2),size(img,1)]);
len = size(img,4);
% �����ȡͼƬ
for iP=1:len
    picture = img(:,:,:,iP); % ����ͼƬ
    picture_fy3_resize = imresize(picture,1); % ��ͼƬ�����ز���
    movie_i(iP) = im2frame(picture_fy3_resize); % ͼƬ����Ϊ������֡
end
% ��fig�в��Ŷ���������1�飬�ٶ�1֡/��
movie(fig,movie_i ,1,5)

%% ������Ƶ��������
 writerObj =VideoWriter(fullfile(savepath,filename)); % ����һ��avi����
 writerObj.FrameRate=5; % ����avi�����Ĳ���������֡����
 open(writerObj); % ��avi����
 writeVideo(writerObj,movie_i); % ������Ķ���д�뵽��Ƶ�ļ���
 close(writerObj); % �رն���
 close(fig);
%  
%  a=img(:,:,:,30);
%  figure,imshow(a,[min(a(:)) max(a(:))/10])