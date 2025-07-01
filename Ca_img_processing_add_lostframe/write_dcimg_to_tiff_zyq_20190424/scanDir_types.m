function [ files ] = scanDir_types( root_dir,file_str)

files={};
if root_dir(end)~='/'
 root_dir=[root_dir,'/'];
end
if file_str(1)~='.'
 file_str=['.' file_str];
end
fileList=dir([root_dir]);  %扩展名
n=length(fileList);
cntpic=0;
for i=1:n
    if strcmp(fileList(i).name,'.')==1||strcmp(fileList(i).name,'..')==1
        continue;
    else
        %fileList(i).name
        if ~fileList(i).isdir % 如果不是目录则跳过
            [~,~,ext] = fileparts(fileList(i).name);
            if (strcmpi(file_str,ext))
            full_name=[root_dir,fileList(i).name];            
                 cntpic=cntpic+1;
                 files(cntpic,1)={full_name};
%              end
            end
        else
            files=[files;scanDir_types([root_dir,fileList(i).name],file_str)];
        end
    end
end

end