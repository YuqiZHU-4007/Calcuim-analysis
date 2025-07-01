function [fpaths, nframes] = ask_for_images(ext)

%% load images
if nargin < 1
    ext = '.tif';
end

imdir = uigetdir;
listdir = dir([imdir '/*' ext]);

% fnames = {listdir.name};
% fpaths = fnames;
% for ii = 1:length(fnames)
%     fpaths{ii} = [imdir '/' fnames{ii}];
% end

% movind = cellfun(@(x) str2num(x(1:end-4)), fnames) + 1;  % order not right
% movind = sort(movind);
nframes = length(listdir);

fpaths = cellfun(@(x) [imdir '/' x], {listdir.name}, 'UniformOutput', false);

%fpaths = cell(nframes, 1);
%for fi = 1:nframes    % caution: fmt dependent
%    fpaths{fi} = [imdir '/' num2str(fi-1) '.bmp'];
%end

