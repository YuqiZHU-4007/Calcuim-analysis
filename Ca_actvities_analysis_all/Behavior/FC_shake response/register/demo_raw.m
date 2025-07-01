% 2016-03-29
% jiangjian@ion.ac.cn
clear,clc;
tic
io.regPath='F:\DUlab\FC_analyse\Ca_actvities_analyze\Behavior\FC_shake response\register\imageJ.jar';
io.ijPath='F:\DUlab\FC_analyse\Ca_actvities_analyze\Behavior\FC_shake response\register\ij.jar';
javaaddpath(io.regPath);
javaaddpath(io.ijPath);
a=imageJ.align();

source='I:\1.Test US\2.Tail free¡ª¡ªData from 117\20210101\fish1\Substack.tif'; %Ä£ï¿½ï¿½ï¿½Ä¼ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Í¼Æ¬ï¿½ï¿½Ã¿ï¿½ï¿½ï¿½Ë¿ï¿½ï¿½Ô¸ï¿½ï¿½ï¿½ï¿½Ô¼ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½É¡ï¿?
targetFolder='I:\1.Test US\2.Tail free¡ª¡ªData from 117\20210101\fish1\substacks/'; %ï¿½Ä¼ï¿½ï¿½Ð£ï¿½Ã¿Ö¡Í¼Æ¬ï¿½ï¿½ï¿½æµ¥ï¿½ï¿½ï¿½Ä¼ï¿½ï¿½ï¿½ï¿½ï¿½jpg/tiffï¿½È¸ï¿½Ê½
outputName='I:\1.Test US\2.Tail free¡ª¡ªData from 117\20210101\fish1\substacks-reg/';
mkdir(outputName)

%ï¿½ï¿½ï¿½ß³Ì·ï¿½ï¿½ï¿½
a.alignMultiImage(source,targetFolder,outputName);

 %ï¿½ï¿½ï¿½ß³Ì·ï¿½ï¿½ï¿½
threadNum=10;
a.alignMultiImageParallel(source,targetFolder,outputName,threadNum);
a.getAlignPoints(source,source)
toc
