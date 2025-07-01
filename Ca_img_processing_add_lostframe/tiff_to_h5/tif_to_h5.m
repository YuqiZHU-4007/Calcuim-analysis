Dir =uigetdir('G:\','tifpath');%'K:\ÐÐÎªtest\20180425_fish2\z02'; % insert path to tiff stack here
h5Dir=uigetdir(Dir,'h5path');
listdir=dir([Dir '/*.*']);
tifpath={};kk=1;
for ii=1:length(listdir)
    if strcmp(listdir(ii).name,'.')==1 || strcmp(listdir(ii).name,'..')==1
        continue;
    else
        if ~isFolder([listdir(ii).folder '\' listdir(ii).name])
            continue;
        else
            tifpath{kk}=[listdir(ii).folder '\' listdir(ii).name];kk=kk+1;
        end
    end
end
tic;
for ii=1:length(tifpath)
    nam=tifpath{ii}    
    h5path=[h5Dir '\z' num2str(ii,'%02d') '.hdf5'];
%     if ~exist([h5Dir '\z' num2str(ii,'%2d')], 'dir') %isFolder([h5Dir '\z' num2str(ii,'%2d')])
%         mkdir([h5Dir '\z' num2str(ii,'%2d')])
%     end
    datasetname='/tifset';
    %Y = read_file(nam);
    Y = seqread(nam);
    %listdir=dir([nam '/*.tif']);
    % Y=imread([nam '\' listdir(1).name]);
    % datasize=size(Y);
    %Y=uint16(zeros(datasize(1),datasize(2),length(listdir)));
    % for ii=1:length(listdir)
    %     Y(:,:,ii)=imread([nam '\' listdir(ii).name]);
    % end
    Datatype=class(Y);
    datasize=size(Y);
    chunksize=[datasize(1) datasize(2) 1];
    if exist(h5path)
        delete(h5path);
    end
    
    h5create(h5path,datasetname,datasize,'Datatype',Datatype,...
        'chunksize',chunksize);
    h5disp(h5path)
    h5write(h5path,datasetname,Y);
    h5disp(h5path)
     
    %     Y2 = read_file(h5path);
    %     for j = 1:datasize(3)
    %         com(:,j)=isequal(Y(:,:,j), Y2(:,:,j));
    %     end
    
    % for ii=1:datasize(3)
    % h5create(h5path,[ '/' num2str(ii)],datasize(1:2),'Datatype',Datatype);
    % h5write(h5path, ['/' num2str(ii)], Y(:,:,ii));
    % end
    %
    % for j = 1:datasize(3)
    %       data = Y(:,:,ii);
    %       start = [1 1 j];
    %       count = [2048 2048 1];
    %       h5write('F:\DUlab\FC analyse\CaImAn-MATLAB-master\²âÊÔ\test.h5','/tifdata',data,start,count);
    % end
end
toc;
