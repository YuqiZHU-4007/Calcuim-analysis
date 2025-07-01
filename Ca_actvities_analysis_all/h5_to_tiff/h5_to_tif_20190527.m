file='I:\3.Juvenile reference brain\¾ÉÄ£°å\Huc_Vmat\20210226fish3_20dpf\20210226fish3_ch2_morph_vol';
nam=[file,'.h5'];
h5disp(nam)
%Y = seqread(nam);
Y = h5read(nam,'/vol');
tic;
zpath = checkpath([file]);
seqwrite(Y, zpath, 'tif')
toc;

