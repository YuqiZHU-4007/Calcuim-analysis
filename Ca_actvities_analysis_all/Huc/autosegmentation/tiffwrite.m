function tiffwrite(imtiff, impath)

for ii = 1:size(imtiff, 3)
    %imwrite(tiff_stack(:,:,:,ii) , 'new_stack.tif' , 'WriteMode' , 'append') ;
    imwrite(imtiff(:,:,ii), impath, 'WriteMode' , 'append');
end