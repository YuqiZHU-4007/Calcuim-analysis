function imio = ImageIO(basepath, T, varargin)
% imio = ImageIO(basepath, T, t0 = 1, varargin)
% 
% use structure instead of class to keep compatibility
% have following methods:
%   readframe(imio, ti, zi)
%   readvolume(imio, ti)
%   readslice(imio, zi)
%   writezt(imio, writepath)
%   writetz(imio, writepath)
% 
% basepath can be these filetypes:
%   'dcimg'  .dcimg file
%   'tif'    folder with all the image files
%   'zt'     folder with subfolders for different z
%   'vol'    folder with multi-page tif for different t
% 
% do not support infer parameters from the filenames so far,
% the filenames should be dealed to the same length
%
% detect image files solely by extentions (for efficient consideration),
% suppose the base image type is tif (tiff) 
% so make sure the image files are actually image files
% 

ip = inputParser;
ip.addRequired('basepath');
ip.addRequired('T');
ip.addOptional('t0', 1);
ip.addParamValue('iscache', false);
ip.addParamValue('nframes', 0);
ip.parse(basepath, T, varargin{:});
inputs = ip.Results;

% if nargin < 5
%     iscache = false;
% end
% 
% if nargin < 3
%     t0 = 1;
% end
% 
% if nargin < 2
%     T = 45;
% end

%% init
imio.basepath  = basepath;
imio.datatype = '4d';
imio.filetype = ''; 
imio.T = T;
imio.t0 = inputs.t0;  t0 = imio.t0;
imio.nframes  = inputs.nframes;
imio.nt       = 0;
imio.nz       = 0;
% imio.t        = [];
% imio.z        = [];
imio.readpath = {};  % ti, zi
imio.iscache  = inputs.iscache;
imio.cache    = {};

%% infer filetype from the path
[pathstr, name, ext] = fileparts(basepath);
if strcmp(ext, '.dcimg') && exist(basepath, 'file')  % dcimg
    [frame, nframes] = dcimgmatlab(0, basepath);
    if imio.nframes == 0  % special treat for the incorrect nframes read from dcimgmatlab
        imio.nframes = nframes;
    end
    nframes = imio.nframes;
    imio.nt = floor(double(nframes - t0 + 1) / T);  % caution: may depend on the format
    imio.nz = T;                                    % nz is not equal to depth
    imio.filetype = 'dcimg';
    
elseif any(strcmp(ext, {'.tif', '.tiff'})) && exist(basepath, 'file')  % single tif, perhaps 2d
    nframes = length(imfinfo(basepath));
    imio.nframes = nframes;
    imio.nt = floor(double(nframes - t0 + 1) / T);  % caution: may depend on the format
    imio.nz = T;
    imio.filetype = 'mptiff';
    
elseif exist(basepath, 'dir')  				       % folder
    
    listdir = dir(basepath);
    listdir = listdir(3:end);  % caution: any better way?
    fnames = {listdir.name};
    fpaths = cellfun(@(x) [basepath '/' x], fnames, 'uniformoutput', false);
    isdir  = [listdir.isdir];
    isimg = cellfun(@(x) strcmp(x(max(1,end-3):end), '.tif') | strcmp(x(max(1,end-4):end), '.tiff'), fnames);
%     ftypes = cellfun(@imftype, fnames, 'uniformoutput', false);
%     isimg  = ~cellfun(@isempty, ftypes);
%     imfmt = ftypes{find(isimg, 1)};
    
    if any(isdir) && any(isimg)     % folder contains both images and subfolders
        error('unable to infer the file type of the input path: folder contains both images and subfolders')

    elseif any(isimg)  % image files
        
%         if ~all(cellfun(@(x) strcmp(x, imfmt), ftypes(isimg)))
%             error('images with different formats in the folder');
%         end

        % plain tif or vol mptif or tif with distinguishable names
        imgnames = fnames(isimg);
        iminfo = imfinfo(fpaths{find(isimg, 1)});
        if length(iminfo) == 1     % plain tif
            imio.filetype = 'tif';
            imio.nframes = sum(isimg);
            imio.nz = T;
            imio.nt = floor(double(imio.nframes - t0 + 1) / T);
            imio.readpath = imgnames;
            
        elseif length(iminfo) > 1  % multi-page
            imio.filetype = 'vol';
            imio.nt = sum(isimg);
            imio.nz = length(iminfo);
            imio.nframes = imio.nz * imio.nt;
            imio.readpath = imgnames;
        end
        
    elseif any(isdir)  % subfolders
       % each z in the folders
       imio.filetype = 'zt';
       imio.nz = sum(isdir);
       %if nz ~= T - 1, error(''); end
       dirnames = fnames(isdir);
       dirpaths = fpaths(isdir);
       for zi = 1:imio.nz
           %imgnames = dir(dirnames{zi});
           %imgnames = imgnames(3:end);  % caution: any better way?
           imgnames = dir([dirpaths{zi} '/*.tif*']);
           imgnames = {imgnames.name};
           imgnames = cellfun(@(x) [dirnames{zi} '/' x], imgnames, 'uniformoutput', false);
           %isimg    = find(~cellfun(@(x) isempty(imftype(x)), imgnames));
           isimg = 1:length(imgnames);
           if zi == 1
               imio.nt = length(isimg);
               imio.readpath = cell(imio.nt, imio.nz);
           elseif length(isimg) ~= imio.nt
               error('subfolders contain different number of images');
           end
           imio.readpath(:,zi) = imgnames(isimg);
       end
       imio.nframes = imio.nz * imio.nt;
       
    else
        error('unable to infer the file type of the input path');
    end
    
else
    error('invalid input path');
end

end

