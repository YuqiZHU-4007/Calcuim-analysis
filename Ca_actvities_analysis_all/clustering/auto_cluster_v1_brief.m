M_0=[];%n*t
M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
%% clustering
numK=[];
ind=index_all;%find(index_all==2);
ind=ind([1:20:length(ind)]);
[gIXx,~,numK,~]=find_best_k_in_range(M_0(ind,:),3:15);
length(unique(gIXx))
cIXx=ind;cIX_reg = (1:size(M_0,1))';%ind;
if isempty(numK)
    numK=20;
end
masterthres=0.6;
para=struct;
para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
for ii=1
[cIX,gIX] = AutoClustering(cIXx,gIXx,M_0,cIX_reg,1,para,1,masterthres);
end
length(unique(gIX))
result_cluster_raw.cIX=cIX;result_cluster_raw.gIX=gIX;
result_cluster_raw.index_all=index_all;
%% some plot
clrmap_name = 'hsv_new';
clrmap = GetColormap(clrmap_name,max(gIX));
length(unique(gIX))
vi='on';
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
outputpath=checkpath([savepath_all '\figures_masterthres' num2str(masterthres) '_' num2str(iii,'%02d')]);
[h,~, numU] = hierplot_zyq_20190530(cIX_iii,gIX_iii,M_0(cIX_iii,:));
B=M_0;
[h]=pushbutton_popupplot_Callback(B,cIX_iii,gIX_iii,clrmap_iii,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
