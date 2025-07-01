%% 0 to 0 CR 
ind=find(CR_ind_up(:,1)==1);ind=intersect(ind,CS_related_hab_ind);
ind2=find(CR_ind_down(:,1)==1);ind2=intersect(ind2, CS_related_hab_down_ind);
indb=union(ind,ind2);indb=setdiff(1:size(A,2),indb);

ind=intersect(find(CR_ind_up(:,2)==1),CS_related_tst_ind);%ind=setdiff(ind,CSmotor_related_tst_ind );
ind2=intersect(find(CR_ind_down(:,2)==1), CS_related_test_down_ind);
ind3=union(ind,ind2);ind3=setdiff(1:size(A,2),ind3);
inter=intersect(indb,ind3);
figure,
scatter3(supervoxel(inter,1),supervoxel(inter,2),supervoxel(inter,3),10,'r','filled');
axis equal;legend('0-1 CR ')
figure,
a=[0 0.05];
patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
    [min(a) min(a) max(a) max(a)],...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
plot(squeeze(mean(A_r(:,[trial.hab(2):trial.hab(2)+5],inter),3)))%trial.test(2):trial.test(3)
figure,imagesc(A(:,inter)',[0 0.05])
%clustering
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
M=A_r(frame_ind,:,inter);M=reshape(M,length(frame_ind)*trial.total,[])';
M_norm=normalize(M,2,'zscore');M=M_norm;
%cluster
numK=[];[gIXx,~,numK,~]=find_best_k_in_range(M,15:40);
cIXx=1:length(inter);cIX_reg = (1:size(M,1))';
masterthres=0.6;
para=struct;para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',20);
for ii=1
[cIX_0to0,gIX_0to0] = AutoClustering(cIXx,gIXx,M,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX_0to0))
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX_0to0));
vi='on';
for iii=1
    cIX_iii=cIX_0to0;gIX_iii=gIX_0to0;clrmap_iii=clrmap;
    %     [h,gIX2, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M)
    %     length(unique(gIX2))
    %2
    B=M;
    [h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
    a=[];
    for ii=1:length(unique(gIX_0to0))
        a(:,ii)=mean(A(:,inter(cIX_0to0(find(gIX_0to0==ii)))),2);
        length(find(gIX_0to0==ii))
    end
    figure,sepplot(1:size(A,1),a,GetColormap('hsv_new',max(gIX_0to0)));hold on;
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
    
    figure,
    indd_0to0={};aa=[];indd_0to0_regressor={};
    for ii=1:length(unique(gIX_0to0))
        M_all=A_r(frame_ind,:,:);M_all=reshape(M_all,length(frame_ind)*trial.total,[]);
        [stimcorr,~] = MotorSourceCorrelation(M_all',mean(M(cIX_0to0(find(gIX_0to0==ii)),:),1),[]);
%         if ii==3 | ii==4 | ii==5 | ii==12 | ii==13
%             n=4;
%         else
%             n=3;
%         end
       n=3;
        thr=mean(stimcorr)+n*std(stimcorr)
        %figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
        indd_0to0{ii}=find(stimcorr>thr);
        indd_0to0_regressor{ii}=inter(cIX_0to0(find(gIX_0to0==ii)));
        aa(:,ii)=mean(A(:,indd_0to0{ii}),2);
        %figure,plot(a(:,ii));hold on;plot(mean(A(:,indd{ii}),2),'r');legend('regressor','stim output')
        scatter3(supervoxel(indd_0to0{ii},1),supervoxel(indd_0to0{ii},2),supervoxel(indd_0to0{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    figure,
    for ii=1:length(unique(gIX_0to0))
        scatter3(supervoxel(indd_0to0_regressor{ii},1),supervoxel(indd_0to0_regressor{ii},2),supervoxel(indd_0to0_regressor{ii},3),10,clrmap(ii,:,:),'filled');hold on;axis equal;
    end
    figure,sepplot(1:size(A,1),aa,GetColormap('hsv_new',max(gIX_0to0)));hold on;
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
    
end


