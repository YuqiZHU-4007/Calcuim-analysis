function h=DrawTiledPics_zyq_20190530(cIX,gIX,absIX,CellXYZ,anat_stack,clrmap)
% isRefAnat = getappdata(hfig,'isRefAnat');
% if ~isRefAnat,

    radius_xy = 16;%10%7;

% initialize image stack
anat_stack2 = zeros([size(anat_stack),3]);
nPlanes = size(anat_stack,3);
dimv_yxz = size(anat_stack);

% make cell-shaped circular mask
% radius_xy = 7;
circlemaskIX = MakeCircularMask(radius_xy,dimv_yxz(1:2));


% clrmap = hsv(round(double(numK)*1.1));
% set transparency
alpha = ones(size(cIX))*0.4;

%% main: coloring of stack
cIX_abs = absIX(cIX);
M_xyz = CellXYZ(cIX_abs,:);
stack_alpha = 0.25;
for i = 1:nPlanes,
    anat_plane = stack_alpha*repmat(squeeze(imNormalize99(anat_stack(:,:,i))),[1 1 3]);
    IX = find(M_xyz(:,3)==i);
    if ~isempty(IX),
        anat_stack2(:,:,i,:) = DrawMasksInRGB(anat_plane,M_xyz(IX,[1,2]),circlemaskIX,clrmap,gIX(IX),alpha(IX));
    end
end

%% 
% if isShowMasks,  
%     nMasks = length(Msk_IDs);
%     % set params
%     cmap2 = jet(nMasks);
%     msk_alpha = 0.3;
%     white_alpha = 0.1;
%     
%     for i = 1:nMasks,
%         msk = full(MASKs.MaskDatabase(:,Msk_IDs(i)));
%         mask_3D = reshape(msk, [MASKs.height, MASKs.width, MASKs.Zs]);
%         
%         % Y-X
% %         masks_XY = max(mask_3D,[],3);
% %         masks_fit = logical(imresize(masks_XY,[size(anat_yx,1),size(anat_yx,2)]));
%         ixs = find(masks_fit);
%         anat_stack2(ixs) = (cmap2(i,1)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % R;
%         ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2);
%         anat_stack2(ixs) = (cmap2(i,2)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % G;
%         ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2)*2;
%         anat_stack2(ixs) = (cmap2(i,3)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % B;
%     end
% end
%% view stack sequentially
% figure;
% for i_plane = 1:nPlanes,
%     im = squeeze(ave_stack2(:,:,i_plane,:));
%     image(im);
%     axis image; axis off
%     pause(0.1)
% end

% lay out

h = figure('Position',[10,20,round(1136*1.4),round(683*1.4)],'color','k');
hold on;

axes('Position',[0,0,1,1]); % BS fix to get the figure background to stay when saving pics
image(zeros(3,3,3));axis off

numRow = 3;
numP = nPlanes - mod(nPlanes,numRow);
numCol = numP/numRow;

for i_plane = 1:numP,
    im = squeeze(anat_stack2(:,:,i_plane,:));
    im = imresize(im,0.25);
    
    [col, row] = ind2sub([numCol,numRow],i_plane); % this is intensionally inverted...
    posVec = [(col-1)*1/numCol,1-row*1/numRow,1/numCol,1/numRow];
    axes('Position',posVec);
    image(im);view([-90 90]);
    axis image; axis off
    %     pause(0.1)
end

end