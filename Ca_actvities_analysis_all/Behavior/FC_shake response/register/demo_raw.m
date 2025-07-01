% 2016-03-29
% jiangjian@ion.ac.cn
clear,clc;
tic
io.regPath='F:\DUlab\FC_analyse\Ca_actvities_analyze\Behavior\FC_shake response\register\imageJ.jar';
io.ijPath='F:\DUlab\FC_analyse\Ca_actvities_analyze\Behavior\FC_shake response\register\ij.jar';
javaaddpath(io.regPath);
javaaddpath(io.ijPath);
a=imageJ.align();

source='I:\1.Test US\2.Tail free����Data from 117\20210101\fish1\Substack.tif'; %ģ���ļ�������ͼƬ��ÿ���˿��Ը����Լ�����������ɡ�?
targetFolder='I:\1.Test US\2.Tail free����Data from 117\20210101\fish1\substacks/'; %�ļ��У�ÿ֡ͼƬ���浥���ļ�����jpg/tiff�ȸ�ʽ
outputName='I:\1.Test US\2.Tail free����Data from 117\20210101\fish1\substacks-reg/';
mkdir(outputName)

%���̷߳���
a.alignMultiImage(source,targetFolder,outputName);

 %���̷߳���
threadNum=10;
a.alignMultiImageParallel(source,targetFolder,outputName,threadNum);
a.getAlignPoints(source,source)
toc
