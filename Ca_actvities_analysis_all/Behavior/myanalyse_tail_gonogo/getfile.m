function getfile
disp('Pick all file expect "settings_constants"');
[filename, pathname] = uigetfile( {'*.mat'}, 'Pick all file expect "settings_constants"','MultiSelect','on');
if size(filename,2)~=0
    for i = 1: length(filename)
        load(strcat(pathname,filename{i}));
    end
else
    load([pathname filename]);
end