%compare 
a1=load('E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\merge_cluster\clust_masterthres_0.6.mat');
a2=load('E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\results_20190731mat');
ind_mean_block=find(a1.result.gIX==21 | a1.result.gIX==29);
ind_1_block=find(a2.result.gIX_rep==11);
index_mean=[];
for ii=unique(index_all)
    ind=find(index_all==ii);
    index_mean=[index_mean;index_all(ind)' [1:length(ind)]'];
end
%index_1=a;
loc1=index_mean(a1.result.cIX(ind_mean_block),:);
loc2=index_1(IX_pass,:);a2.result.cIX_rep(ind_1_block);

ind_this_type_ol=[];
ind_this_type_1=[];ind_this_type_2=[];
ind_this_type_ol_index_all1=[];ind_this_type_1_index_all1=[];ind_this_type_2_index_all1=[];
ind_this_type_ol_index_all2=[];ind_this_type_1_index_all2=[];ind_this_type_2_index_all2=[];
for ii=unique(loc1(:,1))'
    ind1=find(loc1(:,1)==ii);
    ind2=find(loc2(:,1)==ii);
    [C,ia,ib]=intersect(loc1(ind1,2),loc2(ind2,2));
    ind_this_type_ol=[ind_this_type_ol;loc1(ind1(ia),:)];loc2(ind2(ib),:);
    ind_this_type_ol_index_all1=[ind_this_type_ol_index_all1;a1.result.cIX(ind_mean_block(ind1(ia)))];
    ind_this_type_ol_index_all2=[ind_this_type_ol_index_all2;IX_pass(ind2(ib))];
    
    ind_this_type_1=[ind_this_type_1;loc1(ind1(setdiff(1:length(ind1),ia)),:)];
    ind_this_type_1_index_all1=[ ind_this_type_1_index_all1;a1.result.cIX(ind_mean_block(ind1(setdiff(1:length(ind1),ia))))];
    
    ind_this_type_2=[ind_this_type_2;loc2(setdiff(ind2,ind2(ib)),:)];
    ind_this_type_2_index_all2=[ ind_this_type_2_index_all2;IX_pass(setdiff(ind2,ind2(ib)))];
end

ind_this_type_11=loc1;
for ii=1:size(ind_this_type_ol,1)
    ind=find(ind_this_type_11(:,1)==ind_this_type_ol(ii,1) & ind_this_type_11(:,2)==ind_this_type_ol(ii,2));
    ind_this_type_11(ind,:)=[];
end
ind_this_type_22=loc2;
for ii=1:size(ind_this_type_ol,1)
    ind=find(ind_this_type_22(:,1)==ind_this_type_ol(ii,1) & ind_this_type_22(:,2)==ind_this_type_ol(ii,2));
    ind_this_type_22(ind,:)=[];
end
%A=ind_this_type_22;B=ind_this_type_2;a=2;
%setdiff(A(:,a),intersect(B(:,a),A(:,a)))
C=[];
for ii=1:size(ind_this_type_22,1)
    ind=find(index_mean(:,1)==ind_this_type_22(ii,1) & index_mean(:,2)==ind_this_type_22(ii,2));
    %correct
    %C(ii)=sum(index_mean(ind,:)==ind_this_type_22(ii,:));
    if ~isempty(ind)
        ind_this_type_2_index_all1(ii,1)=ind;
%     else
%         ind_this_type_2_index_all1(ii)=[];
    end
end
for ii=1:size(ind_this_type_11,1)
    ind=find(index_1(:,1)==ind_this_type_11(ii,1) & index_1(:,2)==ind_this_type_11(ii,2));
    %correct
    %C(ii)=sum(index_1(ind,:)==ind_this_type_11(ii,:));
    if ~isempty(ind)
        ind_this_type_1_index_all2(ii,1)=ind;
%     else
%         ind_this_type_1_index_all2(ii)=ind;
    end
end
[N_ol,center_ol]= hist( ind_this_type_ol(:,1),1:9);
[N_1,center_1] = hist( ind_this_type_1(:,1),1:9);
[N_2,center_2] = hist( ind_this_type_2(:,1),1:9);
[N,center]= hist( index_mean(:,1),1:9);
%N_ol=N_ol/length(ind_this_type_ol_index_all);
%N_1= N_1/length(ind_this_type_1_index_all);
%N_2= N_2/length(ind_this_type_2_index_all);
figure,bar(1:9,[N_ol;N_1;N_2]'./repmat(N',1,3)); 
xlim([1 8]);
legend('Overlap','Clustered by avg. per block in acq','Clustered by 1.trial per block in acq,hab,test')
%[~,inter1,inter2]=intersect(loc1(:,2),loc2(:,2),'rows');loc1(inter1,:);
clrmap=[0.2 0.2 1];trace=struct;
k = get(0, 'screensize');
colorCS=[0.5 0.5 0.9];
area_win_acq=frame.cs_start:frame.us_start-1;%frame.us_start
type_cal_area='event';%type={'lowest','lowest_all','polyfit','first','mean','sum','event'};
ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;%取取时间段的baseline作为参考判断event
area_win_hab_tst=frame.cs_start:frame.cs_end-1;
for hh=1:3
    switch hh
        case 1
            trace=a2.result.trace_all;
%             trace.hab=reshape(act_all_hab(:,1,:),frame.per_cycle,[]);%cat(2,act_all_hab(:,1,:),act_all_acq_block_mean,act_all_tst(:,1,:));
%             trace.acq=act_all_acq;%reshape(act_all_acq_block_mean,frame.per_cycle*trial.acq_block_num,[]);
%             trace.test=reshape(act_all_tst(:,1,:),frame.per_cycle,[]);
            loc=[ind_this_type_ol_index_all2;ind_this_type_1_index_all2];
            %loc=a1.result.cIX(ind_mean_block);
        case 2
            trace=a2.result.trace_all;
            loc=IX_pass;
            %loc=a2.result.cIX_rep(ind_1_block);
        case 3
            trace=a2.result.trace_all;
            loc=[ind_this_type_ol_index_all2];
    end
    figure,set(gcf,'position',[50,200,2500,400])
    a=[-0.05 0.1];trace_seed=[];
    for ii=1:3
        switch ii
            case 1
                subplot(1,7,1),x=trace.hab(:,loc);isus=false;%set(gca,'position',[50 50 0.1*k(3) 0.85*k(4)]);
            case 2
                subplot(1,7,2:6),x=trace.acq(:,loc);isus=true;%set(gca,'position',[50 50 0.8*k(3) 0.85*k(4)]);
            case 3
                subplot(1,7,7),x=trace.test(:,loc);isus=false;%set(gca,'position',[50 50 0.1*k(3) 0.85*k(4)]);
        end
        patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_start:frame.per_cycle:size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        %plot(x(:,cIX_iii(ind)),'linewidth',0.5,'color',[0.5 0.5 0.5]);hold on
        if isus
            l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
        end
        plot(mean(x,2),'linewidth',2,'color',clrmap);hold on
        ylim(a);xlim([1 size(x,1)])
        box on;
        trace_seed=[trace_seed;mean(x,2)];
    end
    
    area_mean=struct;area_sum_mean=struct;
    for ii=1:3
        switch ii
            case 1
                x=a2.result.trace_all.hab(:,loc);%x=mean(x,2);figure,plot(x)
                area_mean.hab=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
                area_sum_mean.hab(:,1)=mean(area_mean.hab);area_sum_mean.hab(:,2)=std(area_mean.hab)/sqrt(length(loc));
            case 2
                x=a2.result.trace_all_block.acq(:,loc);
                y=reshape(x,frame.per_cycle,trial.acq_block_trial,trial.acq_block_num,[]);%y(:,:,:,1);
                y=mean(y,2);y=reshape(y,frame.per_cycle,trial.acq_block_num,[]);
                for j=1:5 %size(y,3);
                    xx=reshape(y(:,j,:),frame.per_cycle,[]);%xx=mean(xx,2);figure,plot(xx)
                    area_mean.acq{j,1}=calculate_integrate_dfdf(xx,area_win_acq,type_cal_area,ref_win);
                    area_sum_mean.acq(j,1)=mean(area_mean.acq{j,1});area_sum_mean.acq(j,2)=std(area_mean.acq{j,1})/sqrt(length(loc));
                end
            case 3
                x=a2.result.trace_all.test(:,loc);
                area_mean.test=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
                area_sum_mean.test(:,1)=mean(area_mean.test);area_sum_mean.test(:,2)=std(area_mean.test)/sqrt(length(loc));
        end
    end
    area_loc{hh,1}=area_mean;area_loc{hh,2}=area_sum_mean;
end
for hh=1:3
[~,p] = ttest(area_loc{hh,1}.acq{1},area_loc{hh,1}.acq{2})
end

clrmap=[[0.9,0.2,0.9];[0.3,0.8,0];[0,0.3,0.9]];%绿，红,蓝
figure,
scatter([1],[1],10,clrmap(1,:));hold on
scatter([2],[2],10,clrmap(2,:));hold on
scatter([3],[3],10,clrmap(3,:));hold on
for ii=1:1
    switch ii
        case 1
            loc_1=index_1([ind_this_type_1_index_all2],:);
            loc_2=index_1([ind_this_type_2_index_all2],:);
            loc_ol=index_1([ind_this_type_ol_index_all2],:);
    end
    for kk=unique([loc_1(:,1);loc_2(:,1);loc_ol(:,1)])'
        ind_event_in_this_fish_1=loc_1(find(loc_1(:,1)==kk),:);num_1=size(ind_event_in_this_fish_1,1);
        ind_event_in_this_fish_2=loc_2(find(loc_2(:,1)==kk),:);num_2=size(ind_event_in_this_fish_2,1);
        ind_event_in_this_fish_ol=loc_ol(find(loc_ol(:,1)==kk),:);num_ol=size(ind_event_in_this_fish_ol,1);
        load(envbatch{kk})
        num=[1*ones(num_1,1);2*ones(num_2,1);3*ones(num_ol,1)];
        ind_event_in_this_fish=[ind_event_in_this_fish_1;ind_event_in_this_fish_2;ind_event_in_this_fish_ol];
        h=DrawTiledPics_zyq_20190530(ind_event_in_this_fish(:,2),num,[1:size(env.supervoxel,1)],[env.supervoxel(:,2) env.supervoxel(:,1) env.supervoxel(:,3)],env.vol,clrmap);
    end
end