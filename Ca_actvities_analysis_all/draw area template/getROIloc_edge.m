function [MASK,ind_roi]=getROIloc_edge()
[name,path]=uigetfile('/home/ZYQ/W/data/fear conditioning_ZYQ/Fear conditioning/20210101/fish1/*.txt','txt');
[Iname,Ipath]=uigetfile([path '/*.tif'],'.tif');
img=rgb2gray(imread([Ipath Iname]));

[roiname,~,roix,roiy,roiz] = textread([path name],'%s%s%f%f%f','headerlines',1);
n=unique(roiname);
X=1:1:size(img,2);Y=1:1:size(img,1);
[X,Y] = meshgrid(X,Y);X=reshape(X,1,[]);Y=reshape(Y,1,[]);
ind_roi={};
MASK=zeros(size(img));
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
    MASK(ind)=1;
end
figure,imshow(MASK)

 