%% mapback
im=double(env.vol(:,:,1:end-1));rescalegd(im2double(env.vol), [1/10000 1/10000]);%showspv=zeros(size(im,1),size(im,2),3,size(im,3));
showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
%showspv=zeros(size(im,1),size(im,2),3,size(im,3));
% for zz=1:size(showspv,3)
%     showspv(:,:,:,zz)=insertShape(im(:,:,zz), 'filledcircle', [0, 0, 0],'color','black');
% end
indd=ind3;[indd_1toX,indd_0to1];clr=GetColormap('hsv_new',size(indd,2));kk=1;
for ii=1:length(indd)
    id=indd{ii};
    zi=unique(env.supervoxel(id,3));        
    for zz=zi'
      pti =  id(find(env.supervoxel(id,3)==zz));
      slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1],'color',clr(kk,:));
      showspv(:,:,:,zz)=slice;
    end
    kk=kk+1;
end
seqwrite(showspv, 'F:\DUlab\Progress report\202103\0_to_1');

aa=[];idd=[];idd_g=[];
for ii=1:18
    idd=[idd,indd{ii}];
    idd_g=[idd_g;ii*ones(length(indd{ii}),1)];
    aa(:,ii)=mean(A(:,indd{ii}),2);
end
figure,sepplot(1:size(A,1),aa,clr);hold on;
y=[-12 2];
patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))],...
    [min(y) max(y)],'color','r','linestyle','--','linewidth',1);hold on;
xlim([1 size(A,1)])
frame_ind=frame.cs_start:frame.us_start-1;
M=A_r(frame_ind,:,:);M=reshape(M,length(frame_ind)*trial.total,[])';
M_norm=normalize(M,2,'zscore');M=M_norm;
[h]=pushbutton_popupplot_Callback(M,idd,idd_g,clr,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
