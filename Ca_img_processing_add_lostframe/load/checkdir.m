function dirpath = checkdir(dirpath)
% dirpath = checkdir(dirpath)

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end