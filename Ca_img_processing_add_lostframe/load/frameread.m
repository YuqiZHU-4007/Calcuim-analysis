function [frame, nframes] = frameread(imgpath, fid)
% [frame, nframes] = frameread(imgpath, fid)

[frame, nframes] = frmread(imgpath, fid(1));

end

function [frame, nframes] = frmread(imgpath, fid)

[pathstr, name, ext] = fileparts(imgpath);

if strcmp(ext, '.dcimg')                              % dcimg
    [frame, nframes] = dcimgmatlab(fid-1, imgpath);
    frame = frame';
    clear functions;
elseif exist(imgpath, 'dir')  						  % folder
    fnames = dir([imgpath '/*.tif']);  % assume it's tif
    fnames = {fnames.name};
    nframes = length(fnames);
    if length(unique(cellfun(@length, fnames))) == 1  % dealed names
        for fi = 1:nframes
            fnames{fi} = [imgpath '/' fnames{fi}];
        end
        frame = imread(fnames{fid});
    else  											  % pure number names
        frame = imread([imgpath '/' num2str(fid) '.tif']);
    end
else  												  % assume it's prefix name of tif
	frame = imread([imgpath num2str(fid) '.tif']);
    nframes = length(dir(pathstr)) - 2;
end

end

