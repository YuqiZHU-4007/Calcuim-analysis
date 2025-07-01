clc;clear all
[name,path]=uigetfile('I:\1.Test US\2.Tail free¡ª¡ªData from 117\20210101\fish1\*.txt','txt');
[Iname,Ipath]=uigetfile([path '\*.tif'],'.tif');
img=rgb2gray(imread([Ipath Iname]));

[roiname,~,roix,roiy,roiz] = textread([path name],'%s%s%f%f%f','headerlines',1);
n=unique(roiname);
X=1:1:size(img,2);Y=1:1:size(img,1);
[X,Y] = meshgrid(X,Y);X=reshape(X,1,[]);Y=reshape(Y,1,[]);
ind_roi={};
for ii=1:length(n)
    ind=find(strcmp(string(roiname),string(n(ii))));
    x=roix(ind);
    y=roiy(ind);
    
    mask=insertShape(img,'Polygon',{reshape([x,y]',1,[])},'linewidth',3);
    [in,on]=inpolygon(X,Y,x,y);numel(X(in))
    numel(X(on))
    ind=sub2ind([size(img,1),size(img,2)],Y(in),X(in));
    mask2=img;mask2(ind)=0;
    
    figure,subplot(1,3,1),imshow(img);hold on;
    subplot(1,3,2),imshow(mask);
    subplot(1,3,3),imshow(mask2);
    ind_roi{ii}=ind;
end


% xv =ROI2.t_No(ind);
% yv = ROI2.XY(ind);
% xq = 1:1:size(img,2);
% yq = 1:1:size(img,1);
% [xq,yq] = meshgrid(xq,yq);xq=reshape(xq,1,[]);yq=reshape(yq,1,[]);
% 
% [in,on] = inpolygon(xq,yq,xv,yv);numel(xq(in))
% figure
% plot(xv,yv) % polygon
% hold on
% plot(xq(~in),yq(~in),'bo') 
% plot(xq(in),yq(in),'r+') % points inside
% %% points outside
% hold off
% xlim([1 640]);ylim([1 512]);
