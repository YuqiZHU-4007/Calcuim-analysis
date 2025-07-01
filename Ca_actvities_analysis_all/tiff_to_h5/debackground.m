
%env_path='C:\Users\pc\Desktop\standard brain\light-sheet\env.mat';%'K:\行为test\20180425_fish2\env.mat';
% [envname envpath]=uigetfile('E:\A_Data_lightsheet\Data_huc\*.mat','*.mat');
% load([envpath envname]);
[volname volpath]=uigetfile([envpath '*.tiff'],'*.tiff');
%nam='C:\Users\pc\Desktop\standard brain\light-sheet\20190606_fish2_huc.tiff';%'K:\行为test\20180425_fish2\Plane11_stack.hdf5';
%h5disp(nam);
%Y = seqread(nam);
Y = read_file([volpath volname]);
dim=size(Y);
layernum=11;
Datatype=class(Y);
%subtract mean of each line out of mask
Y_de=Y;
for ii=1:dim(3)
    %mask=env.volmask(:,:,ii);
    %mask = imdilate(mask, strel('disk', env.width/16+1));%%膨胀mask
    %figure,imshow(mask);
    %     ind=find(~mask);
    %     sub=ind2sub([dim(1) dim(2)],ind);
    frame=Y(:,:,ii);
    %figure,imshow(frame,[min(frame(:)) max(frame(:))]);
    %     frame(:,1000:end)=0;figure,imshow(frame,[min(frame(:)) max(frame(:))]);
    %     mask(:,1000:end)=0;figure,imshow(mask);
    %frame(mask==1)=0;
    m=sum(frame,2)./sum(frame~=0,2);
    %m=mean(frame,2,'omitnan');
    frame=Y(:,:,ii)-uint16(repmat(m,1,dim(1)));
    %a(ii,1)=mean(frame(~mask));
    %aa(ii,1)=mean(frame(mask));
    %     frame(~mask)=1;
%     figure,
%     subplot(1,2,1),imshow(Y(:,:,ii),[min(min(Y(:,:,ii))) max(max(Y(:,:,ii)))]);
%     subplot(1,2,2),imshow(frame,[min(frame(:)) max(frame(:))]);

    Y_de(:,:,ii)=frame;
end
[~,name,format]=fileparts(volname);
tiffwrite(Y_de, [volpath name '_Deback_line.tiff']);

% %subtract mean of area out of mask
% Y_de=Y;
% for ii=1:dim(3)
%     mask=env.volmask(:,:,ii);
%     %mask = imdilate(mask, strel('disk', env.width/16+1));%%膨胀mask
%     %figure,imshow(mask);
%     %     ind=find(~mask);
%     %     sub=ind2sub([dim(1) dim(2)],ind);
%     frame=Y(:,:,ii);
%     %figure,imshow(frame,[min(frame(:)) max(frame(:))]);
%     %     frame(:,1000:end)=0;figure,imshow(frame,[min(frame(:)) max(frame(:))]);
%     %     mask(:,1000:end)=0;figure,imshow(mask);
%     m=mean(frame(~mask));
%     frame=Y(:,:,ii)-uint16(m);
%     b(ii,1)=mean(frame(~mask));
%     bb(ii,1)=mean(frame(mask));
%     %     frame(~mask)=1;
% %     figure,
% %     subplot(1,2,1),imshow(Y(:,:,ii),[min(min(Y(:,:,ii))) max(max(Y(:,:,ii)))]);
% %     subplot(1,2,2),imshow(frame,[min(frame(:)) max(frame(:))]);
%     Y_de(:,:,ii)=frame;
% end
%tiffwrite(Y_de, 'C:\Users\pc\Desktop\standard brain\light-sheet\20190606_fish2_huc_Deback.tiff');

% h5path='K:\行为test\20180425_fish2\Plane11_stack_debackground.hdf5';
% datasetname='/tifset';
% datasize=dim;
% chunksize=[dim(1) dim(2) 1];
% h5create(h5path,datasetname,datasize,'Datatype',Datatype,...
%     'chunksize',chunksize);
% h5disp(h5path)
% h5write(h5path,datasetname,Y_de);
% h5disp(h5path)

