clear all;
Dir =uigetdir('G:\','tifpath');% insert path to tiff stack here
% '/mnt/X/data/fear_conditioning/huc/20190514/fish3/rec25106988/raw tiff'
mkdir([Dir,'_h5path']);
h5Dir=[Dir,'_h5path'];
listdir=dir([Dir '/z01/*.tif']);
y= imread([listdir(1).folder,'/',listdir(1).name]);datasize=size(y);Datatype=class(y);
T=length(listdir);
l=dir([Dir '/*z*']);N_plane=length(l)-1;
datasetname='/default';chunksize=[datasize(1) datasize(2) 1];
compress=4;

tic;
for ti = 1:T;
    Y=zeros(datasize(1),datasize(2),N_plane);Y=uint16(Y);
    for zi = 1:N_plane
        yy=imread([Dir '/z',num2str(zi,'%02d'),'/',num2str(ti,'%04d'), '.tif']);
        Y(:,:,zi)=yy;
    end
    h5path=fullfile([h5Dir '/t' num2str(ti,'%04d'),'.h5']);
    if exist(h5path)
        delete(h5path);
    end
    h5create(h5path,datasetname,[datasize(1) datasize(2) N_plane],'Datatype',Datatype,...
        'chunksize',chunksize,'Deflate',compress);
    h5write(h5path,datasetname,Y);
    toc;
end
toc;

%% compare
ti=200;
filename=fullfile([h5Dir '/t' num2str(ti,'%04d'),'.h5']);
data_h5 = h5read(filename,datasetname);
%show_spv_GUI(data_h5);
data_tif=zeros(datasize(1),datasize(2),N_plane);data_tif=uint16(data_tif);
   for zi = 1:N_plane
        yy=imread([Dir '/z',num2str(zi,'%02d'),'/',num2str(ti,'%04d'), '.tif']);
        data_tif(:,:,zi)=yy;
   end
% show_spv_GUI(data_tif);   
length(find((data_h5-data_tif)~=0))
%% 
%show_spv_GUI(Y);
% poolobj = parpool(poolsize);
% spmd(poolsize)
%     ti_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize T]);    
%     for ii = 1:length(ti_in_this_thread)
%         ti = ti_in_this_thread(ii);
%         [labindex ti]     
%         Y=repmat(y,1,1,N_plane);
%         for zi = 1:N_plane
%             yy=imread([Dir '/z',num2str(zi,'%02d'),'/',num2str(ti,'%04d'), '.tif']);
%             Y(:,:,zi)=yy;
%         end
%         h5path=fullfile([h5Dir '/t' num2str(ti,'%04d'),'.h5']);
%         if exist(h5path)
%             delete(h5path);
%         end
%         h5create(h5path,datasetname,[datasize(1) datasize(2) N_plane],'Datatype',Datatype,...
%             'chunksize',chunksize,'Deflate',compress); 
%         h5write(h5path,datasetname,Y);
%         toc;
%     end
% end
% delete(gcp('nocreate'));


