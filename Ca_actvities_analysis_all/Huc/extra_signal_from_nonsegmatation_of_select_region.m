function [ROI_traces,para] = extra_signal_from_nonsegmatation_of_select_region
clc;clear all;close all;
% automatically chose ROIs
% 201411_YangChen
%20141201_cyin
%20160831_cyin: for autoROI H2B-G6f nucleus expressing data; only differ in detect_by_ring_2d_n.m subfunction
%20190102_zyq

%%import vol
[fn_tif,fn_path] = uigetfile('G:/*.tiff','Vol'); % using tif stack file
%%import registration image
image_dir = uigetdir('G:/','registration image dir'); % using tif stack file
%read stack
stack = my_tiffread([fn_path,fn_tif]);%imread(fn_tif); %my_tiffread.m same as imread.m for tiff stack frame by frame
%outpath
outpath=checkpath([fn_path,'subregoin_supervol']);
%%  parameters
poolsize=3;
blocksize = ceil(2550/ poolsize);
para.loc_nucleus=true;%logic
para.isregistration=true;
para.minrad = 3;
para.maxrad = 8;
para.thres = 0;
radii =  para.minrad: para.maxrad;

para.region_num=input('input the number of regions you selected to extra signal\n');
region_name=cell(para.region_num,1);
outpath_region=cell(para.region_num,1);
edge_x=cell(para.region_num,1);edge_y=cell(para.region_num,1);
mask_region=cell(para.region_num,1);
for zz=1:para.region_num
    region_name{zz}=input('input the name of regions you selected to extra signal\n','s');
    para.T{zz}=input(['input the layer of ', region_name{zz},' regions, eg. 2:5\n']);
    outpath_region{zz}=checkpath([outpath '\' region_name{zz}]);
    %manually define ROI area
    for ii= para.T{zz}
        slice =stack(:,:,ii);
        slice = normim(slice);
        gs = fspecial('gaussian', opt.maxrad, 1);
        slice = imfilter(slice, gs);
        
        h1=figure('name',['layer' num2str(ii) '-' region_name{zz}]); imshow(slice);
        [BW,c,r]= roipoly(slice);
        if ~isempty(BW)
            mask_region{zz}(:,:,ii) = BW; % manually define ROI area
            edge_x{zz}{ii,1}=c;edge_y{zz}{ii,1}=r;
        else
            mask_region{zz}(:,:,ii)=zeros(size(slice));
            edge_x{zz}{ii,1}=0;edge_y{zz}{ii,1}=0;
        end
        h2=figure; imshow(mask_region{zz}(:,:,ii));
        close(h1);close(h2);
    end
end

%segmentation

for zz=1:para.region_num
f=cell(length(para.T),1);
supervoxel=[];
for ii= para.T
    slice =stack(:,:,ii);
    slice = normim(slice);
    %slice_back = slice;
    % smooth  smooth width should be smaller than cell diameter
    gs = fspecial('gaussian',  para.maxrad*2+1, 1);%
    slice = imfilter(slice, gs);
    %slice = double(slice);
    
    if para.loc_nucleus
        %[points, fim, radmap] = detect_by_ring_2d_n(slice, radii);
        [points, fim, radmap] = detect_by_ring_2d_n(slice, radii,para.thres, mask_region{zz}(:,:,ii));
    else
        [points, fim, radmap] = detect_by_ring_2d_cytosol(slice, radii, para.thres, mask_region{zz}(:,:,ii));
    end
    
    [y, x] = find(points);
    rad = sum(radmap, 3) / 2;
    rad = rad(points);
    
    rin = radmap(:,:,1);
    rin = rin(points);
    rex = radmap(:,:,2);
    rex = rex(points);
    scores = fim(points);
    if ~isempty(x)
        idx = remove_overlap([x y], rin, scores);
        rad = rad(idx);
        rin=rin(idx);
        rex = rex(idx);
        x = x(idx);
        y = y(idx);
    end
    
    ind = sub2ind(size(slice), y, x);%将下标转换为线性索引
    points = false(size(points));
    points(ind) = true;
    
    tshow = slice;
    tshow(points) = true;
    % show = slice / max(slice(:));
    % tshow = show;
    % tshow(points) = true;
    
    f{ii,1} = figure('Name','autoROI','color', 'white', 'position', [400   200   size(slice)]);% 1 plot
    imshow(tshow);hold on;
    %text(x, y, num2str((1:length(x))'),'color','y');
    text(x, y, num2str((1:length(x))'),'fontsize',11,'HorizontalAlignment','center');hold on;
    plot(edge_x{zz}{ii,1},edge_y{zz}{ii,1},'r','linewidth',1.5);hold on;
    ellipse(rad,rad,zeros(size(rad)),x,y,1.5);
    saveas(f{ii,1},[outpath,'\',num2str(ii),'.tif']);
    saveas(f{ii,1},[outpath,'\',num2str(ii),'.fig']);
    close(f{ii,1});
    supervoxel = [supervoxel; [x y] ii*ones(length(rin),1) rin rex];
end
para.supervoxel{zz,1}=supervoxel;
end

%%%% time series of each ROI
%supervol region
spvzind = cell(length(para.T),1);
for zi = para.T
    spvzind{zi} = find(supervoxel(:,3) == zi)';
end
nspv = size(supervoxel, 1);
spvregion = cell(nspv, 1);
for pi = 1:nspv
    xp = supervoxel(pi,1);
    yp = supervoxel(pi,2);
    rp = floor(supervoxel(pi,4));  % use outer radius?
    [x, y] = meshgrid(xp-rp-1:xp+rp+1,yp-rp-1:yp+rp+1);
    ind = find((x - xp).^2 + (y - yp).^2 <= rp.^2);
    r = y(ind); c = x(ind);
    try
        ind = sub2ind(size(slice), r, c);%sub2ind函数将矩阵中指定元素的行列下标转换成存储的序号，即线性索引号。
    catch
        r(r<1) = 1; r(r>height) = height;
        c(c<1) = 1; c(c>width) = width;
        ind = sub2ind(size(slice), r, c);
        ind = unique(ind);
    end
    spvregion{pi} = ind;
end
para.spvregion=spvregion;
para.spvzind=spvzind;
%extral signal
ROI_traces=cell(length(para.T),1);
if length(para.T)>=2
    for ii= para.T
        filepath=[image_dir,'/z',num2str(ii, '%.2d')];
        listdir = dir([filepath '/*.tif']);
        poolobj = parpool(poolsize);
        spmd(poolsize)
            ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize length(listdir)]);
            actth = zeros(length(spvzind{ii}), length(ti_in_this_thread));
            for jj = 1:length(ti_in_this_thread)
                ti = ti_in_this_thread(jj);
                [labindex ti]
                filename = listdir(ti).name;
                image_registration=imread([filepath,'\',filename]);
                for ppi = spvzind{ii}
                    actth(ppi,jj) = mean(image_registration(spvregion{ppi}));
                end
            end
        end
        ROI_traces{ii,1}= cat(2, actth{:});
        delete(poolobj);
    end
else
    for ii= para.T
        filepath=[image_dir,'/z',num2str(ii, '%.2d')];
        listdir = dir([filepath '/*.tif']);
        %     filename = listdir(1).name;
        %     slice=imread([filepath,filename, '.tif']);
        %     image_registration=[];
        %     image_registration=zeors(size(slice,1),size(slice,2),length(listdir));
        for jj = 1:length(listdir)
            filename = listdir(jj).name;
            image_registration=imread([filepath,filename, '.tif']);
            for ppi = spvzind{ii}
                ROI_traces{ii,1}(ppi, jj) = mean(image_registration(spvregion{ppi}));
            end
        end
    end
end
save([fn_path '/ROI_traces.mat'], 'ROI_traces', 'para', '-v7.3');

end



