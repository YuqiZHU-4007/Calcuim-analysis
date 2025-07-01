function seqwrite(img, imgpath, fmt)
% seqwrite(img, imgpath, fmt = 'tif')

if nargin < 3
    fmt = 'tif';
end

checkpath(imgpath);

if length(size(img)) == 3
    for ii = 1:size(img, 3)
        try
            imwrite(img(:,:,ii), [imgpath '/' num2str(ii, '%04d') '.' fmt]);
        catch 
            tiffwrite(img(:,:,ii), [imgpath '/' num2str(ii, '%04d') '.' fmt]);
        end
    end
elseif length(size(img)) == 4 && size(img, 3) == 3        % assume it's rgb
    for ii = 1:size(img, 4)
        imwrite(img(:,:,:,ii), [imgpath '/' num2str(ii, '%04d') '.' fmt]);
    end
elseif length(size(img)) == 4 && size(img, 4) == 3        % will be wrong if depth is 3
    warning('will be wrong if depth is 3');
    for ii = 1:size(img, 3)
        imwrite(squeeze(img(:,:,ii,:)), [imgpath '/' num2str(ii, '%04d') '.' fmt]);
    end
elseif length(size(img)) == 2
    imwrite(img, [imgpath '/0001.' fmt]);
else
    error('invalid input');
end

