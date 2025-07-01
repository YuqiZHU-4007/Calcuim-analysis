Dir =uigetdir('G:\','tifpath');%'K:\ÐÐÎªtest\20180425_fish2\z02'; % insert path to tiff stack here
h5Dir=uigetdir(Dir,'h5path');
listdir=dir([Dir '/*.*']);
tifpath={};kk=1;
for ii=1:length(listdir)
    if strcmp(listdir(ii).name,'.')==1 || strcmp(listdir(ii).name,'..')==1
        continue;
    else
        if ~isfolder([listdir(ii).folder '\' listdir(ii).name])
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
    %     if ~exist([h5Dir '\z' num2str(ii,'%2d')], 'dir') %isfolder([h5Dir '\z' num2str(ii,'%2d')])
    %         mkdir([h5Dir '\z' num2str(ii,'%2d')])
    %     end
    datasetname='/tifset';
    %Y = read_file(nam);
    Y = seqread(nam);
    Datatype=class(Y);
    datasize=size(Y);
    batch_width=10;
    batch=100;
    batch=max(min(round(batch/batch_width)*batch_width,datasize(3)),1);
    chunksize=[datasize(1) datasize(2) batch];
    if exist(h5path)
        delete(h5path);
    end
    
    h5create(h5path,datasetname,datasize,'Datatype',Datatype,...
        'chunksize',chunksize);
    h5disp(h5path)
    h5write(h5path,datasetname,Y);
    h5disp(h5path)
end
toc;

%listdir=dir([nam '/*.tif']);
% Y=imread([nam '\' listdir(1).name]);
% datasize=size(Y);
%Y=uint16(zeros(datasize(1),datasize(2),length(listdir)));
% for ii=1:length(listdir)
%     Y(:,:,ii)=imread([nam '\' listdir(ii).name]);
% end