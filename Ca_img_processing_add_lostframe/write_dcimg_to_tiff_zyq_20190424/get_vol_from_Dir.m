function get_vol_from_Dir
T=25;
file_str='tiff';
file_name='vol';
fliedir=uigetdir('path of vol folder');
outpath=uigetdir('path of save vol folder');

imgpathes=scanDir_types( fliedir,file_str);
imgpathes_vol={};kk=1;
for ii=1:length(imgpathes)
    p=imgpathes{ii};
    if strcmp(p(end-length(file_str)-3:end-length(file_str)-1),file_name)
        imgpathes_vol{kk,1}=p;kk=kk+1;
    else
        continue;
    end
end
nbatch = length(imgpathes_vol);
stack=imread(imgpathes_vol{1},T);
dim=size(stack);
all_vol_stack_T=uint16(zeros(dim(1),dim(2),nbatch));
for ii=1:nbatch 
    p=imgpathes_vol{ii};
%     copyfile(p,fullfile(outpath,['vol_' num2str(ii,'%04d') '.' file_str]))
    stack=imread(p,T);
    all_vol_stack_T(:,:,ii)=stack;
end

seqwrite(all_vol_stack_T,fullfile(outpath,'all_vol_stack_T'),'tif');
save([outpath '\vol_path'],'imgpathes_vol');