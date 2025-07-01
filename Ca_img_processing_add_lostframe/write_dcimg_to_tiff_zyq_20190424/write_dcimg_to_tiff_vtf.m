setenv('MW_MINGW64_LOC','Z:\Ca_img_processing_add_lostframe\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup

T=25;
fmt='tif';
file_str='dcimg';

fliedir=uigetdir('path of dcimg folder');
imgpathes=scanDir_types(fliedir,file_str);
nbatch = length(imgpathes);
%elapsedTime=zeros(nbatchi,1);
for batchi=1:nbatch
    imgpath=imgpathes{batchi};
    %[imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
    savepath=imgpath(1:end-6);%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\';
    
    zpath={};
    for jj=1:T
        zpath{jj,1} = checkpath([savepath '/raw tiff/z' num2str(jj, '%02d')]);
    end
    imio = dcimgreader(imgpath, T, 1);  %, opt.t0); % auto estimate t
    t=imio.nt;
    for ii=1:t
        frame=imio.readvolume(ii);
        for  jj=1:T
            imwrite(frame(:,:,jj), [zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt]);
            %     frame_lzw=imread([zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt]);
            %     com(ii,jj)=isequal(frame(:,:,jj), frame_lzw);
        end
    end
end