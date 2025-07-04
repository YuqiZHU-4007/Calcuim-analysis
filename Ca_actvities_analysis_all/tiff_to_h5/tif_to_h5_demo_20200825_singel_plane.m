clear all;
Dir =uigetdir('G:\','tifpath');% insert path to tiff stack here
% '/mnt/X/data/fear_conditioning/huc/20190514/fish3/rec25106988/raw tiff'
[~,filename]=fileparts(Dir)
mkdir([Dir,'_h5']);
h5Dir=[Dir,'_h5'];
listdir=dir([Dir '/*.tif']);
y= imread([listdir(1).folder,'/',listdir(1).name]);datasize=size(y);Datatype=class(y);
T=length(listdir);
datasetname='/default';chunksize=[datasize(1),datasize(2) 1];
compress=4;

h5path=fullfile([h5Dir '.h5']);
if exist(h5path)
    delete(h5path);
end

h5create(h5path,datasetname,[datasize(1),datasize(2) T],'Datatype',Datatype,...
    'chunksize',chunksize,'Deflate',compress);
h5disp(h5path); 

for ii = 1:T
  1;
  frame=imread(fullfile(Dir ,'/',[filename,'_',num2str(ii,'%01d'), '.tif']));
  h5write(h5path,datasetname,frame,[1 1 ii],[datasize(1),datasize(2),1]);
end
