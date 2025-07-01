file1=fopen('/data/Plane11.stack');
A=fread(file1,'uint16');
A_r=reshape(A,2048,1188,[]);
a=A_r(:,:,11);
figure,imshow(a,[min(a(:)) max(a(:))])
h5path='/data/plane11.hdf5'
datasetname='/plane11';
datasize=size(A_r);
chunksize=[datasize(1) datasize(2) 1];
Datatype='uint16';
h5create(h5path,datasetname,datasize,'Datatype',Datatype,...
    'chunksize',chunksize);
h5disp(h5path)
h5write(h5path,datasetname,A_r);
h5disp(h5path)