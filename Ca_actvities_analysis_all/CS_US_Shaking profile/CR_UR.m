colorCS=[0.5 0.5 0.9];
%% CR in hab.
A=stimcorr.CS(:,7);%A=A(find(index_all==6));
B=act_all_hab;
B=act_all_tst;%B=B(:,:,find(index_all==6));
isus=false;
figure,hist(A,100);
thr=mean(motorcorr.sep_spon(:,1),'omitnan')+2*std(motorcorr.sep_spon(:,1),[],1,'omitnan');
figure,hist(motorcorr.sep_spon(:,1),100);
p_hab_base=nan(length(index_all),1);h_hab_base=nan(length(index_all),1);
for kk=1:length(index_all)
    y=mean(B(frame.cs_start:frame.cs_end,:,kk),1);
    x=mean(B(1:frame.cs_start-1,:,kk),1);
    [p_hab_base(kk,1),h_hab_base(kk,1) ]= ranksum(x,y);
    %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
end
[~,q] = mafdr(p_hab_base);h_hab_base_adj=nan(size(q));h_hab_base_adj(find(q<0.05),1)=1;
disp(length(find(h_hab_base==1)));
disp(length(find(h_hab_base_adj==1)));

ind_CR={};ind_CR_cutmotor={};
for zz=1
    figure;a=[-0.005 0.2];kk=1;b=0.5;[-0.2:0.1:0.4];
    for tt=b
        %ind=find(A>tt);
        ind=find(sum(h_hab_base_adj,2)>0.8*size(h_hab_base_adj,2));%ind=find(A>tt & sum(h_hab_base_adj,2)>0.8*size(h_hab_base_adj,2));
        switch zz
            case 1
                ind_CR{kk}=ind;
            case 2
                disp(['cut motor ' num2str(length(find(motorcorr.sep_spon(ind,1)>thr)))])
                ind(find(motorcorr.sep_spon(ind,1)>thr))=[];
                ind_CR_cutmotor{kk}=ind;
        end 
        subplot(length([b]),1,kk),
        x=B(:,:,ind);x=reshape(x,[],1,length(ind));
        patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_start:frame.per_cycle:size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        if isus
            l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
        end
        plot(mean(x,3),'linewidth',1.5);hold on;
        shadedErrorBar([1:size(x,1)],x,{@(x) mean(x,3)',@(x) std(x,[],3)},'lineprops',{[1 0.2 0.2]},'transparent',1,'patchSaturation',0.1); hold on
        ylim(a);title(num2str(tt));
        kk=kk+1;
        disp(['num ' num2str(length(ind)) ' sd: ' num2str(max(std(x,[],3)))]);
    end
    fraction_in_region_in_clust_CR=[];num_in_fish_CR=[];
    for ii=unique(index_all(ind))'
        ind_in_this_fish=(find(index_all==ii));
        i=find(index_all(ind)==ii);
        load(envbatch{ii});load([savebatch{ii} '\brain arearegion_mask.mat']);
        figure,scatter3(env_all.supervoxel(ind(i),2),env_all.supervoxel(ind(i),1),env_all.supervoxel(ind(i),3),'filled');axis equal;
        xlim([1 env.height]);ylim([1 env.width])
        %showspv  = mapback_show_zyq_20190729(A(ind(i)), env_all.supervoxel(ind(i),:),1:length(i),env.vol,'dot');
        %h=show_spv_GUI(showspv);
        [~,a,~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,ones(size(ind(i))),ind(i),...
            [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
            [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],[1 0.2 0.2],true);
        fraction_in_region_in_clust_CR(:,ii)=a(:,1);
        num_in_fish_CR(ii)=length(i)/length(ind_in_this_fish);
        title(['Fish ' num2str(ii,'%02d')]);
    end
    figure,bar(unique(index_all(ind)),num_in_fish_CR(unique(index_all(ind))));ylim([0 1]);
end
%figure,hist(A(find(index_all==6)),100)
%% UR
ind_UR={};isus=true;
p_acq_base=nan(length(index_all),size(act_all_acq_block_mean,2));h_acq_base=nan(length(index_all),size(act_all_acq_block_mean,2));
for ii=1:size(act_all_acq_block_mean,2)
    figure,hist(stimcorr.US(:,ii),100);xlim([-1 1]);
    for kk=1:length(index_all)
        y=act_all_acq_block_mean(frame.us_start:frame.cs_end,ii,kk);
        x=act_all_hab_mean(1:frame.cs_start-1,1,kk);
        [p_acq_base(kk,ii),h_acq_base(kk,ii) ]= ranksum(x,y);
    end
end
h_acq_base_adj=nan(length(index_all),size(act_all_acq_block_mean,2));
for ii=1:size(act_all_acq_block_mean,2)
    [~,q] = mafdr(p_acq_base(:,ii));h_acq_base_adj(find(q<0.05),ii)=1;
    disp(length(find(h_acq_base(:,ii)==1)));
    disp(length(find(h_acq_base_adj(:,ii)==1)));
end

fraction_in_region_in_clust_UR={};num_in_fish_UR={};
for ss=1:size(act_all_acq_block_mean,2)
    %figure,hist(stimcorr.US(:,ss),100);
    figure;a=[-0.005 0.2];kk=1;b=0.5;[0:0.1:0.6];
    for tt=b
        %ind=find(stimcorr.US(:,ss)<tt);
        ind=find(h_acq_base(:,ss)==1);
        ind_UR{kk,ss}=ind;
        subplot(length([b]),1,kk),
        x=act_all_acq_block_mean(:,ss,ind);x=reshape(x,[],1,length(ind));
        patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_end:frame.per_cycle:size(x,1)]'...
            [frame.cs_start:frame.per_cycle:size(x,1)]']',...
            repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        if isus
            l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
        end
        plot(mean(x,3),'linewidth',1.5);hold on;
        shadedErrorBar([1:size(x,1)],x,{@(x) mean(x,3)',@(x) std(x,[],3)},'lineprops',{[1 0.2 0.2]},'transparent',1,'patchSaturation',0.1); hold on
        ylim(a);%title(num2str(tt));
        title(num2str(length(ind)));
        kk=kk+1;
        disp(['num ' num2str(length(ind)) ' sd: ' num2str(mean(std(x,[],3)))]);
    end
    for ii=unique(index_all(ind))'
        ind_in_this_fish=(find(index_all==ii));
        i=find(index_all(ind)==ii);
        load(envbatch{ii});load([savebatch{ii} '\brain arearegion_mask.mat']);
        %         figure,scatter3(env_all.supervoxel(ind(i),2),env_all.supervoxel(ind(i),1),env_all.supervoxel(ind(i),3),'filled');axis equal;
        %         xlim([1 env.height]);ylim([1 env.width])
        [~,a,~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,ones(size(ind(i))),ind(i),...
            [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
            [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],[1 0.2 0.2],true);
        fraction_in_region_in_clust_UR{ss}(:,ii)=a(:,1);
        num_in_fish_UR{ss}(ii)=length(i)/length(ind_in_this_fish);
        title(['Fish ' num2str(ii,'%02d') 'Acq. ' num2str(ss,'%02d')]);
    end
    figure,bar(unique(index_all(ind)),num_in_fish_UR{ss}(unique(index_all(ind))));ylim([0 1]);
end
%% CS & US convergence
ind=ind_CR{1};num_ind_CS_US={};ind_CS_US={};
clrmap=[1 1 0.2;1 0.2 0.2;0.2 1 0.2];
figure;a=[-0.005 0.1];
for zz=1:3
    subplot(3,1,zz),
    switch zz
        case 1
            iii=intersect(ind_CR{1},ind_UR{1,ss});
        case 2
            iii=setdiff(ind_CR{1},intersect(ind_CR{1},ind_UR{1,ss}));
        case 3
            iii=setdiff(ind_UR{1,ss},intersect(ind_CR{1},ind_UR{1,ss}));
    end
    %x=act_all_acq_block_mean(:,1,iii);x=reshape(x,[],1,length(iii));
    x=act_all_hab_mean(:,1,iii);x=reshape(x,[],1,length(iii));
    patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_start:frame.per_cycle:size(x,1)]']',...
        repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
    plot(mean(x,3),'linewidth',1.5);hold on;
    shadedErrorBar([1:size(x,1)],x,{@(x) mean(x,3)',@(x) std(x,[],3)},'lineprops',{clrmap(zz,:)},'transparent',1,'patchSaturation',0.1); hold on
end
for ii=unique(index_all(ind))'
    ind_in_this_fish=(find(index_all==ii));
    i=find(index_all(ind)==ii);
    load(envbatch{ii});load([savebatch{ii} '\brain arearegion_mask.mat']);
    for ss=1:size(act_all_acq_block_mean,2)
        %          figure,
        %          scatter3(env_all.supervoxel(ind(i),2),env_all.supervoxel(ind(i),1),env_all.supervoxel(ind(i),3),'r','filled');hold on;%CS
        ind2=ind_UR{1,ss};
        i2=find(index_all(ind2)==ii);
        %          scatter3(env_all.supervoxel(ind2(i2),2),env_all.supervoxel(ind2(i2),1),env_all.supervoxel(ind2(i2),3),'g','filled');hold on;%US
        iii=intersect(ind(i),ind2(i2));
        %          scatter3(env_all.supervoxel(iii,2),env_all.supervoxel(iii,1),env_all.supervoxel(iii,3),'y','filled');hold off;axis equal;%CS&US
        num_ind_CS_US{ii,ss}(:,1)=length(iii);
        num_ind_CS_US{ii,ss}(:,2)=length(setdiff(ind(i),iii));
        num_ind_CS_US{ii,ss}(:,3)=length(setdiff(ind2(i2),iii));
        ind_CS_US{ii,ss}{1}=iii;
        ind_CS_US{ii,ss}{2}=setdiff(ind(i),iii);
        ind_CS_US{ii,ss}{3}=setdiff(ind2(i2),iii);
        [~,a,~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,[ones(num_ind_CS_US{ii,ss}(:,1),1);2*ones(num_ind_CS_US{ii,ss}(:,2),1);3*ones(num_ind_CS_US{ii,ss}(:,3),1)],...
            [ind_CS_US{ii,ss}{1};ind_CS_US{ii,ss}{2};ind_CS_US{ii,ss}{3}],...
            [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
            [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap,true);
        title(['Fish ' num2str(ii,'%02d') 'Acq. ' num2str(ss,'%02d')]);
    end
end