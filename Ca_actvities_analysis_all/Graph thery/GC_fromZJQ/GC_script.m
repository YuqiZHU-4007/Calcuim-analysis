% [GC_nonzero_train1,GCstrength_train1]=Groupar_mls_concatenated(graph.node_mean_act_training(:,451:450+75*6)',75);
% [GC_nonzero_train2,GCstrength_train2]=Groupar_mls_concatenated(graph.node_mean_act_training(:,901:450+75*12)',75);
% [GC_nonzero_train2,GCstrength_train3]=Groupar_mls_concatenated(graph.node_mean_act_training(:,1351:450+75*18)',75);
% [GC_nonzero_train2,GCstrength_train4]=Groupar_mls_concatenated(graph.node_mean_act_training(:,1801:450+75*24)',75);
% [GC_nonzero_train2,GCstrength_train5]=Groupar_mls_concatenated(graph.node_mean_act_training(:,2251:450+75*30)',75);
%
% [GC_nonzero_train,GCstrength_train]=Groupar_mls_concatenated(graph.node_mean_act_training',75);
% [GC_nonzero_pre,GCstrength_pre]=Groupar_mls_concatenated(graph.node_mean_act_training(:,1:450)',75);
% [GC_nonzero_test,GCstrength_test]=Groupar_mls_concatenated(graph.node_mean_act_training(:,2701:end)',75);
%
%
% [GC_nonzero_act1,GCstrength_act1]=Groupar_mls_concatenated(graph.node_mean_act',3103);
% [GC_nonzero_act2,GCstrength_act2]=Groupar_mls_concatenated(graph.node_mean_act2',3095);

[GC_nonzero_train1,GCstrength_train1]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.hab(3)*frame.per_cycle+1:trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle)',frame.per_cycle);
[GC_nonzero_train2,GCstrength_train2]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle+1:trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*2)',frame.per_cycle);
[GC_nonzero_train3,GCstrength_train3]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*2+1:trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*3)',frame.per_cycle);
[GC_nonzero_train4,GCstrength_train4]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*3+1:trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*4)',frame.per_cycle);
[GC_nonzero_train5,GCstrength_train5]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*4+1:trial.acq(3)*frame.per_cycle)',frame.per_cycle);

[GC_nonzero_train,GCstrength_train]=Groupar_mls_concatenated(graph.node_mean_act_training',frame.per_cycle);
[GC_nonzero_pre,GCstrength_pre]=Groupar_mls_concatenated(graph.node_mean_act_training(:,1:trial.hab(3)*frame.per_cycle)',frame.per_cycle);
[GC_nonzero_test,GCstrength_test]=Groupar_mls_concatenated(graph.node_mean_act_training(:,trial.acq(3)*frame.per_cycle+1:end)',frame.per_cycle);


[GC_nonzero_act1,GCstrength_act1]=Groupar_mls_concatenated(graph.node_mean_act',3103);
GC_nonzero_act1_spaced={};GC_strength_act1_spaced={};
figure,
for ii=[0:ceil(2*60/fs.ca):size(graph.node_mean_act,2)-ceil(5*60/fs.ca)]
    ind=[1:ceil(5*60/fs.ca)]+ii;
    [GC_nonzero_act1_spaced{ii},GC_strength_act1_spaced{ii}]=Groupar_mls_concatenated(graph.node_mean_act(:,ind)',size(graph.node_mean_act,2));
    plot();
end

[GC_nonzero_act2,GCstrength_act2]=Groupar_mls_concatenated(graph.node_mean_act2',3095);
%% 两两计算GC
graph_GC_paired_two=struct;
for ss=1:7
    switch ss
        case 1
            frame_ind=1:trial.hab(3)*frame.per_cycle;
        case {2,3,4,5,6}
            frame_ind=trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*(ss-2)+1:trial.hab(3)*frame.per_cycle+trial.acq_block_trial*frame.per_cycle*(ss-1);
        case 7
            frame_ind=trial.acq(3)*frame.per_cycle+1:trial.test(3)*frame.per_cycle;
    end
    GC_nonzero=[];GCstrength=[];pF=[];pWald=[];
    for ii=1:size(graph.nodes,1)
        for jj=1:size(graph.nodes,1)
            if ii==jj
                GC_nonzero(ii,jj)=1;GCstrength(ii,jj)=inf;pF(ii,jj)=nan;pWald(ii,jj)=nan;
            else
                x=graph.node_mean_act_training(ii,frame_ind);
                y=graph.node_mean_act_training(jj,frame_ind);
                [a,b,~,c,d]=Groupar_mls_concatenated([x;y]',frame.per_cycle);
                GC_nonzero(ii,jj)=a(1,2);GCstrength(ii,jj)=b(1,2);pF(ii,jj)=c(1,2);pWald(ii,jj)=d(1,2);
            end
        end
    end
    graph_GC_paired_two=setfield(graph_GC_paired_two,['GC_nonzero',num2str(ss)],GC_nonzero);
    graph_GC_paired_two=setfield(graph_GC_paired_two,['GCstrength',num2str(ss)],GCstrength);
    graph_GC_paired_two=setfield(graph_GC_paired_two,['pF',num2str(ss)],pF);
    graph_GC_paired_two=setfield(graph_GC_paired_two,['pWald',num2str(ss)],pWald);
end
save('K:\4.合作_20191022\code\GC_fromZJQ\graph_GC_paired_two.mat','graph_GC_paired_two');
%% show
% C=GC_nonzero_act1;
% [idx,~,~,~] = kmeans(C',5,'Distance','correlation','Replicates',10);
% [~,idx_ind]=sort(idx);
% %C=corr(input.node_mean_act(idx_ind,:)');
% h=figure('name','C');imagesc(C(idx_ind,idx_ind));colormap('hot');colorbar;
figure('name',['GC_nonzero spon']);
set(gcf,'position',[50,200,1100,200]);
subplot(1,4,1),imagesc(GC_nonzero_act1);title('pre spon');
subplot(1,4,2),imagesc(GC_nonzero_act2);title('aft spon');
subplot(1,4,3),imagesc([(GC_nonzero_act2-GC_nonzero_act1) > 0]);title('increased');
subplot(1,4,4),imagesc([(GC_nonzero_act1-GC_nonzero_act2) > 0]);title('decreased');
a=GC_nonzero_act1;b=GC_nonzero_act2;
a=sum(a,1);b=sum(b,1);
%a=(GC_nonzero_act2-GC_nonzero_act1) > 0;a=sum(a,1);
GC_increased_ind=find(b > a);
[~,GC_increased_ind_max]=max(b-a);
figure,histogram(b-a,[0:1:max(b-a)]);
figure,
scatter3(graph.nodes(:,1),graph.nodes(:,2),graph.nodes(:,3),10,[0.5 0.5 0.5],'filled');hold on;
scatter3(graph.nodes(GC_increased_ind,1),graph.nodes(GC_increased_ind,2),graph.nodes(GC_increased_ind,3),12,'k','filled');hold on;
scatter3(graph.nodes(GC_increased_ind_max,1),graph.nodes(GC_increased_ind_max,2),graph.nodes(GC_increased_ind_max,3),20,'r','filled');hold on;
axis equal;
%为什么act2对角线没有因果关系
figure('name',['GC_nonzero training']);
set(gcf,'position',[50,200,1800,200]);
for ss=1:7
    a=[-1 1];
    b=getfield(graph_GC_paired_two,['GC_nonzero',num2str(ss)]);
    subplot(1,7,ss),imagesc(b,a);title([num2str(ss)],'fontsize',12);
end
subplot(1,7,1),imagesc(GCstrength_pre,a);title('hab');
subplot(1,7,2),imagesc(GCstrength_train1,a);title('acq.1');
subplot(1,7,3),imagesc(GCstrength_train2,a);title('acq.2');
subplot(1,7,4),imagesc(GCstrength_train3,a);title('acq.3');
subplot(1,7,5),imagesc(GCstrength_train4,a);title('acq.4');
subplot(1,7,6),imagesc(GCstrength_train5,a);title('acq.5');
subplot(1,7,7),imagesc(GCstrength_test,a);title('test');
%traing阶段刺激会影响GC吗？需不需要自发？