ȱ�� tiffread/write Ҫȷ��matlab·������
ȱ�� pipeline_registration_single_thread��ԭ�ļ���û�У�
ȱ�� show_spv
ȱ�� readfunc����extract_signal 210����
�Ķ�����
dcimgreader���ĺ�����dcimgreader_zyq( )
ref_volume�ĺ�����ref_volume_zyq
pipeline_registration_ori #no change based raw code
readframe
extract_signal�������֡������frameΪnan����back=0����Ӧ backgroundyΪ0��

20180828
����ļ�I:\Ca_img_processing����֡��\func\in use ������Ϊ��������Ҫ������Ƕ�׺�����
I:\Ca_img_processing����֡��\load\frameread.m+ Image4DReader.m +seqread.m+seqwrite.m+seqwrite_rgb.m

20180915���£�ref_volume_zyq,���β�opt.savepath,��������averge volume��ͼ�����nsample�еĵ�zi�������diff<0.8,�����error�� 
%%cut it in20181002 beacuase all fishes in 20180624's corr computed vy ref_volume was smaller than 0.5

20181012����:processing_main,��disp���Ⱥ�ʱ�䣬�Լ���׼��ʾÿ������ʱ��


�汾HCimg4.2.033
HCimage4.2(make sure your C:\Windows\System32\dcimgapi.dll is 14.3.639.4514) 
and also need C/C++ compiler(setup steps: setenv(��the path of your C/C++ compiler��);mex -setup),
now we used MinGw64.

����֡�棨���Կ��а棩
1��t0=0���ӵ�һ֡��ʼ
2��minrad���ھ�����ֵ��/maxrad���⾶����ֵ�����뾶 ������minrad��maxrad��ֵ��������������ھ���minrad��maxrad�����⾶���ھ�+1��maxrad�����
3��setenv('MW_MINGW64_LOC','C:\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');
mex -setup��ֻҪ��Mingw�����ԣ�ע��C��C++�ı�������Ҫ���ã����ú�ر�����matlab��������
���������calllibָ�����Ͳ�����->���Ĭ��·�����������Z:\Ca_img_processing_add_lostframe���ر�����matlab������matlab��

 setenv('MW_MINGW64_LOC','F:\HCImage\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup
top:setenv('MW_MINGW64_LOC','D:\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup
bottom:setenv('MW_MINGW64_LOC','Z:\Ca_img_processing_add_lostframe\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup

4��supervoxel segmentation�����zi�ڼ��� aaa��ʱthresΪ[]��thresΪ�ж�fim���ָ����������������score�����ֵ��fim<thres���point����Ϊһ��ϸ����
4����ʱ��3�ˣ���registration  55.4267h���ڶ���31h���� extra signal 2.4h


