缺少 tiffread/write 要确保matlab路径中有
缺少 pipeline_registration_single_thread（原文件就没有）
缺少 show_spv
缺少 readfunc（在extract_signal 210处）
改动函数
dcimgreader：改后命名dcimgreader_zyq( )
ref_volume改后命名ref_volume_zyq
pipeline_registration_ori #no change based raw code
readframe
extract_signal（如果少帧，读出frame为nan，则back=0，对应 backgroundy为0）

20180828
添加文件I:\Ca_img_processing加少帧版\func\in use （里面为各步骤主要函数的嵌套函数）
I:\Ca_img_processing加少帧版\load\frameread.m+ Image4DReader.m +seqread.m+seqwrite.m+seqwrite_rgb.m

20180915更新：ref_volume_zyq,加形参opt.savepath,用来存算averge volume的图。如果nsample中的第zi层的最大的diff<0.8,则出现error。 
%%cut it in20181002 beacuase all fishes in 20180624's corr computed vy ref_volume was smaller than 0.5

20181012更新:processing_main,加disp进度和时间，以及配准显示每层所用时间


版本HCimg4.2.033
HCimage4.2(make sure your C:\Windows\System32\dcimgapi.dll is 14.3.639.4514) 
and also need C/C++ compiler(setup steps: setenv(‘the path of your C/C++ compiler’);mex -setup),
now we used MinGw64.

加少帧版（测试可行版）
1、t0=0；从第一帧开始
2、minrad（内径估计值）/maxrad（外径估计值）：半径 。会在minrad：maxrad的值中算出差异最大的内径（minrad：maxrad）和外径（内径+1：maxrad）组合
3、setenv('MW_MINGW64_LOC','C:\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');
mex -setup（只要是Mingw都可以，注意C和C++的编译器都要设置，设置后关闭所有matlab，重启）
（如果出错：calllib指针类型不符合->变成默认路径，重新添加Z:\Ca_img_processing_add_lostframe，关闭所有matlab，重启matlab）

 setenv('MW_MINGW64_LOC','F:\HCImage\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup
top:setenv('MW_MINGW64_LOC','D:\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup
bottom:setenv('MW_MINGW64_LOC','Z:\Ca_img_processing_add_lostframe\mingw-w64\x86_64-4.9.2-posix-seh-rt_v3-rev1\mingw64\');mex -setup

4、supervoxel segmentation输出：zi第几层 aaa此时thres为[]（thres为判断fim（分割出的区域的最大内外score差）的阈值，fim<thres则该point不算为一个细胞）
4、用时（3核）：registration  55.4267h（第二次31h）。 extra signal 2.4h


