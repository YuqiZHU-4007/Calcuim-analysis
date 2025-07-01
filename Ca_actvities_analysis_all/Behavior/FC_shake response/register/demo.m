clear,clc;
tic
io.regPath='I:\register\imageJ.jar';
io.ijPath='I:\Fiji_64bits.app\ij.jar';
javaaddpath(io.regPath);
javaaddpath(io.ijPath);
a=imageJ.align();

%和模板(每个trial第一帧)
sourceFolder=uigetdir('Z:\','source folder');
%targetFolder=uigetdir('Z:\','target folder');
outputName=uigetfile('Z:\','output folder');
%sourceFolder='H:\20180509\fish2\roi-fin-total';
%targetFolder='H:\20180509\fish2\roi-fin-total'; 
%outputName='H:\20180509\fish2\roi_fin_reg_result\'; 
%a.alignMultiImage(source,targetFolder,outputName);
threadNum=20;%线程数
%a.alignMultiImageParallel(source,targetFolder,outputName,threadNum);
framenum=1:72000;
trialnum=60;%
totframe_percycle=1200;

%parpool(2)
% if isempty(gcp('nocreate'))==1
%     pool=parpool('local');
% end
%spmd
for ii=1:trialnum
    sourcenum=(ii-1)*totframe_percycle+1;
    source=[sourceFolder '\' num2str(sourcenum) '.tif'];
    for jj=[1:totframe_percycle]+sourcenum-1
        target=[sourceFolder '\' num2str(jj) '.tif'];
        pp{jj}=a.getAlignPoints(source,target);
        delta_r2(ii,jj-sourcenum+1)=2*pp{jj}(1,2);delta_c2(ii,jj-sourcenum+1)=2*pp{jj}(1,1);
        delta_r(jj)=2*pp{jj}(1,2);delta_c(jj)=2*pp{jj}(1,1);
    end
end
%end
%delete(gcp)
toc
% for jj=1:size(pp,2)
%     ind=fix(jj/totframe_percycle)+min(mod(jj,totframe_percycle),1);
%     if mod(jj,totframe_percycle)~=0
%         delta_r2(ind,mod(jj,totframe_percycle))=2*pp{jj}(1,2);
%         delta_c2(ind,mod(jj,totframe_percycle))=2*pp{jj}(1,1);
%     else
%         delta_r2(ind,totframe_percycle)=2*pp{jj}(1,2);
%         delta_c2(ind,totframe_percycle)=2*pp{jj}(1,1);
%     end
% end
h1=figure;subplot(2,1,1),plot(delta_r);title('r');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_r) max(delta_r)],'color','r');hold on;
txt=num2str([1:trialnum]');
text([1:trialnum]'*totframe_percycle-800,(min(delta_r)+max(delta_r))/2*ones(trialnum,1),txt);
subplot(2,1,2),plot(delta_c);title('c');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_c) max(delta_c)],'color','r');hold on;
text([1:trialnum]'*totframe_percycle-800,(min(delta_c)+max(delta_c))/2*ones(trialnum,1),txt);



%和前一帧(每个trial第一帧和自己)
for ii=1:trialnum
   for jj=[1:totframe_percycle]
   source=[sourceFolder '\' num2str(max(jj-1,1)+(ii-1)*totframe_percycle) '.tif']; 
   target=[sourceFolder '\' num2str(jj+(ii-1)*totframe_percycle) '.tif']; 
    pp_bef1{jj+(ii-1)*totframe_percycle}=a.getAlignPoints(source,target);
    delta_r_bef1(jj+(ii-1)*totframe_percycle)=2*pp_bef1{jj+(ii-1)*totframe_percycle}(1,2);
    delta_c_bef1(jj+(ii-1)*totframe_percycle)=2*pp_bef1{jj+(ii-1)*totframe_percycle}(1,1);
    delta_r_bef12(ii,jj)=2*pp_bef1{jj+(ii-1)*totframe_percycle}(1,2);
    delta_c_bef12(ii,jj)=2*pp_bef1{jj+(ii-1)*totframe_percycle}(1,1);
   end
end
toc
h2=figure;subplot(2,1,1),plot(delta_r_bef1);title('r');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_r_bef1) max(delta_r_bef1)],'color','r');hold on;
txt=num2str([1:trialnum]');
text([1:trialnum]'*totframe_percycle-800,(min(delta_r_bef1)+max(delta_r_bef1))/2*ones(trialnum,1),txt);
subplot(2,1,2),plot(delta_c_bef1);title('c');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_c_bef1) max(delta_c_bef1)],'color','r');hold on;
text([1:trialnum]'*totframe_percycle-800,(min(delta_c_bef1)+max(delta_c_bef1))/2*ones(trialnum,1),txt);


[y_3sd,startpoint_sd]=findpeak_use3sd(delta_r_bef1,10);
[y_t,startpoint_t]=findpeak_usetendency(delta_r_bef1);
%reshape
re_startpoint_sd(:,1)=fix(startpoint_sd/totframe_percycle);
re_startpoint_sd(:,1)=mod(startpoint_sd,totframe_percycle);
save([outputName '\agin_result.mat'], 'delta_r', 'delta_r2','delta_c', 'delta_c2', ...
    'delta_r_bef1','delta_r_bef12','delta_c_bef1','delta_c_bef12','pp','pp_bef1',...
     'y_3sd','startpoint_sd','y_t','startpoint_t','-v7.3');