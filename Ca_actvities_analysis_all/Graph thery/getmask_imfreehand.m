%getmask
function mask=getmask_imfreehand(vol,num_region)
I=max(vol,[],3);
mask=zeros(size(I,1),size(I,2));
boundary={};    
h=figure;imshow(I,[min(I(:)) max(I(:))]);hold on;
for ii=1:num_region
    b=[];BW=false(size(I));
    for jj=1:2
        %imagesc(mask);hold on;
        h_draw = imfreehand;
        if numel(h_draw) ~= 0
            BW = BW | createMask(h_draw);
            B = bwboundaries(BW, 'noholes');
            if numel(B) > 0
                b = [b;B{1}];
            end
        end
        boundary{ii,1}=b;
        mask(BW==1)=ii;
    end
    for zz=1:ii-1
        plot(boundary{zz,1}(:,2),boundary{zz,1}(:,1));hold on;
    end
end
close(h);
%figure,imshow(mask,[min(mask(:)) max(mask(:))]);
% addToolbarExplorationButtons(gcf)
end