function ppath = checkpath(ppath)

if ~exist(ppath, 'dir')
    mkdir(ppath);
end