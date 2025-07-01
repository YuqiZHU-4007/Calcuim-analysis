savepathi='K:\7.年会figure\cluster\during consitioning';fmt='fig';
%% single clust plot
num_clust=[23,21,20,10,4];cIX_iii=cIX;gIX_iii=gIX;
clrmap_iii=GetColormap('jet',length(num_clust));
%clrmap_iii=colormap(jet(length(num_clust)));
yli=[-0.005 0.01];kk=1;pt=[];ht=[];
isplotback=false;
%CS patch
color=[1 1 1];
for jj=num_clust
    a=find(gIX_iii(ind)==jj);a=ind(a);
    for ii=2:6
        if ii==1
            B=squeeze(act_all_hab_mean);b=mean(B(:,cIX_iii(a)),2);linevisible=false;linestyle='--';
        elseif ii==7
            B=squeeze(act_all_tst_mean);b=mean(B(:,cIX_iii(a)),2);linevisible=false;linestyle='-';
        else
            B=reshape(act_all_acq_block_mean(:,ii-1,:),50,[]);b=mean(B(:,cIX_iii(a)),2);linevisible=true;
        end
        h=figure;
%         if linevisible %US line
%             l3=plot([frame.us_start frame.us_start]*fs.ca,[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',2);hold on
%         end
        patch1=patch([frame.cs_start,...
            frame.cs_end,...
            frame.cs_end,...
            frame.cs_start]*fs.ca,...
            [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
            color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
        if isplotback %backgroud
            Ba=reshape(act_all_acq_block_mean(:,2,:),50,[]);
            plot(mean(Ba(:,cIX_iii(a)),2),'linewidth',2,'color',[0.5 0.5 0.5]);hold on
        end
        %trace
        c=7;
        %plot(b,'linewidth',3.5,'color',clrmap_iii(kk,:),'LineStyle','--');hold on;
        plot([frame.cs_start-c:frame.us_start-2]*fs.ca,b(frame.cs_start-c:frame.us_start-2),'linewidth',4.5,'color',clrmap_iii(kk,:),'LineStyle','--');hold on;
        plot([frame.cs_start:frame.us_start-2]*fs.ca,b(frame.cs_start:frame.us_start-2),'linewidth',4.5,'color',clrmap_iii(kk,:));hold on;
        set(gca,'visible','off');
        xlim([frame.cs_start-c frame.us_start-2]*fs.ca);ylim(yli);
        saveas(h,[savepathi '\clusnum' num2str(jj,'%02d') '_session' num2str(ii,'%02d') ],fmt);
        %close(h);
        
        t1=mean(B,2);t1=t1(frame.cs_start:frame.us_start-1);t2=b(frame.cs_start:frame.us_start-1);
        %t1=B(frame.cs_start:frame.us_start-1,:);t2=B(frame.cs_start:frame.us_start-1,cIX_iii(ind));
        [ht(kk,ii),pt(kk,ii),~,~]=ttest2(t1,t2,'Vartype','equal');
        if ii==2
            obj = scalebar;
            obj.XLen = 2;
            obj.YLen = 0.005; %X-Length, 10.
            obj.XUnit = 's';            %X-Unit, 'm'.
            %obj.Position = [3.419354838709680,-0.003248688046647];  %move the whole SCALE position.
            obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
            obj.hTextX.FontSize=16;obj.hTextY.FontSize=16;
            obj.hLineX(1).Color=[1 1 1];obj.hLineX(2).Color=[1 1 1];
            obj.hLineY(1).Color=[1 1 1];obj.hLineY(2).Color=[1 1 1];
            obj.hLineX(1).LineWidth=2;obj.hLineX(2).LineWidth=2;
            obj.hLineY(1).LineWidth=2;obj.hLineY(2).LineWidth=2;
            obj.hTextX.Color=[1 1 1];obj.hTextY.Color=[1 1 1];
        end
    end
    kk=kk+1;
end

%% averaged across num_clus
num_clust=[20,21,7,23,10,4,9,5,11,12,13,19,15,17];
%num_clust=[20,21,7,23,10,4,9];%无到有
num_clust=[11,13,12,5,19,15,17];%有到增加,17
%num_clust=[1,2,3,6,8,14,16,18,22];%不变的
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
yli=[-0.02 0.04];kk=1;pt=[];ht=[];ind=[];
isplotback=false;
for jj=num_clust
    ind=cat(1,ind,find(gIX_iii==jj));
end
% % %%correct by ranksum test
h_base=nan(1,max(unique(cIX_iii)))';p_base=nan(1,max(unique(cIX_iii)))';
for kk=unique(cIX_iii)'
    x=squeeze(mean(act_all_hab(1:frame.cs_start-1,:,kk),1));
    y=squeeze(mean(act_all_hab(frame.cs_start:frame.cs_end,:,kk),1));
    [p_base(kk),h_base(kk) ]= ranksum(x,y,'tail','left');
    %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
end
h_base_adj=zeros(1,max(unique(cIX_iii)))';h_tst_adj=zeros(1,max(unique(cIX_iii)))';
[FDR_base,q_base] = mafdr(p_hab_base);h_base_adj(find(q_base<0.05))=1;
[FDR_base,q_base] = mafdr(p_tst_base);h_tst_adj(find(q_base<0.05))=1;

% h_tst=nan(1,max(unique(cIX_iii)))';p_tst=nan(1,max(unique(cIX_iii)))';h_tst_adj=zeros(1,max(unique(cIX_iii)))';
% for kk=unique(cIX_iii)'
%     x=squeeze(mean(act_all_hab(frame.cs_start:frame.cs_end,:,kk),1));
%     y=squeeze(mean(act_all_tst(frame.cs_start:frame.cs_end,:,kk),1));
%     [p_tst(kk),h_tst(kk) ]= ranksum(x,y,'tail','left');
%     %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
% end
% h_acq=nan(trial.acq_block_num,max(unique(cIX_iii)))';p_acq=nan(trial.acq_block_num,max(unique(cIX_iii)))';h_acq_adj=zeros(trial.acq_block_num,max(unique(cIX_iii)))';
% for kk=unique(cIX_iii)'
%     x=squeeze(mean(act_all_hab(1:frame.cs_start-1,:,kk),1));
%     for jj=1:trial.acq_block_num
%         y=squeeze(mean(act_all_acq_block_mean(frame.cs_start:frame.cs_end,jj,kk),1));
%         [p_acq(jj,kk),h_acq(jj,kk) ]= ranksum(x,y,'tail','left');
%     end
%     %[p(kk),h(kk) ]= ttest(x,y,'tail','left');
% end
% FDR_acq=nan(trial.acq_block_num,max(unique(cIX_iii)))';q_acq=nan(trial.acq_block_num,max(unique(cIX_iii)))';
% for jj=1:trial.acq_block_num
%     [FDR_acq(jj,:),q_acq(jj,:)] = mafdr(p_acq(jj,:));h_acq_adj(jj,find(q_acq(jj,:)<0.05))=1;
% end
rank_ind=find(h_base_adj==0);
rank_ind=find(h_base_adj==1 & h_tst_adj==0);rank_ind=find(h_tst_adj==1);rank_ind=find(h_base_adj==0 & h_tst_adj==1);%i_ind=i_ind(500);
[ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(ind));ind=ind(i_ind);
%cur motor related
figure,hist(motorcorr.sep_CS(cIX_iii(ind),7))
ind(find(motorcorr.sep_CS(cIX_iii(ind),7)>=0.3 ))=[];%& motorcorr.sep_spon(cIX_iii(ind),7)>=0.4

% diff
dif=[];
for ii=unique(index_all)'
    num_l=find(index_all==ii);
    dif(ii,1)=length(intersect(num_l,cIX_iii(ind)))/length(num_l);
end
%dif=[bl_fraction;bn_fraction]';dif=sum(dif(num_clust,:));
%per region
fraction_in_region_in_clust_c=[];ii=1;
for kk=unique(index_all)'
    ind_in_this_fish=(find(index_all==kk));
    cIX_in_this_clus=cIX_iii(ind);%
    [IX,i_ind_in_this_fish,i_cIX]=intersect(ind_in_this_fish' ,cIX_in_this_clus');
    unique(index_all(cIX_iii(ind(i_cIX))));
    unique(gIX_iii(ind(i_cIX)));
    %e=load(envbatch{kk});
    load([savebatch{kk} '\brain arearegion_mask.mat']);
    %figure,scatter( env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2), env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1))
    %     [~,fraction_in_region_in_clust(:,:,ii),~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),[1:length(i_cIX)],...
    %         [env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),3)],clrmap_iii(num_clust,:),true);
    [loc_in_region_in_clust,fraction_in_region_in_clust_c(:,:,ii),loc_in_region_cell,id_in_region_cell]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind(i_cIX)),cIX_iii(ind(i_cIX)),...
        [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
        [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap_iii(num_clust,:),false);
    
    title(['Fish ' num2str(kk,'%02d')]);
    ii=ii+1;
end
x=sum(mean(fraction_in_region_in_clust_c(:,num_clust,1:4),3,'omitnan'),2,'omitnan');y=sum(mean(fraction_in_region_in_clust_c(:,num_clust,5:end),3,'omitnan'),2,'omitnan');
plot_fraction_regions([x,y],Lable,[1 0 0.2;0.2 0 1]);ylim([0 2]);a=[x,y];
plot_fraction_regions([x-y]./[x+y],Lable,[1 0.5 0]);ylim([-1 1]);legend off;
%mapback
for kk=[2,3,4,6,12,13]
    load(envbatch{kk});
    num_l=find(index_all==kk);
    [~,~,i_ind]=intersect(num_l,cIX_iii(ind));cmap=repmat([1 0.2 0.2],length(unique(gIX_iii(ind))),1);
    %figure,scatter3(X(cIX_iii(ind(i_ind)),1),X(cIX_iii(ind(i_ind)),2),X(cIX_iii(ind(i_ind)),3),'filled');axis equal;
    DrawTiledPics_zyq_20190530(cIX_iii(ind(i_ind)),gIX_iii(ind(i_ind)),[1:size(A,3)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,cmap);
    %DrawTiledPics_zyq_20190530(cIX_iii(ind(i_ind)),gIX_iii(ind(i_ind)),[1:size(A,3)],[env_all.supervoxel(:,1) env_all.supervoxel(:,2) env_all.supervoxel(:,3)],permute(env.vol,[2 1 3]) ,cmap);
end
%trace
clrtrace='r';clrmap_iii(kk,:);'r';
%CS patch
color=[1 1 1];
for tt=1
    h=figure;a=[-0.01 0.05];
    patch1=patch([frame.cs_start,...
        frame.cs_end,...
        frame.cs_end,...
        frame.cs_start]*fs.ca,...
        [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
        color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
    %     if isplotback %backgroud
    %         plot(mean(B,2),'linewidth',2,'color',[0.5 0.5 0.5]);hold on
    %     end
    for ii=[1,7]
        if ii==1
            B=squeeze(act_all_hab_mean);b=mean(B(:,cIX_iii(ind)),2);linevisible=false;linestyle='--';
        elseif ii==7
            B=squeeze(act_all_tst_mean);b=mean(B(:,cIX_iii(ind)),2);linevisible=false;linestyle='-';
        else
            B=reshape(act_all_acq_block_mean(:,ii-1,:),50,[]);b=mean(B(:,cIX_iii(ind)),2);linevisible=true;
        end
        if linevisible %US line
            l3=plot([frame.us_start frame.us_start]*fs.ca,[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',2);hold on
        end      
        %trace
        %shadedErrorBar([1:frame.per_cycle]',B(:,cIX_iii(ind))',{@(x) mean(x,1),@(x) std(x,[],1)./sqrt(size(x,1))},'lineprops',{'r'},'transparent',1,'patchSaturation',0.1); hold on
        plot([1:frame.per_cycle]*fs.ca,b,'linewidth',3,'color',clrtrace,'linestyle',linestyle);hold on;
        plot([frame.cs_start:frame.us_start]*fs.ca,b(frame.cs_start:frame.us_start),'linewidth',3.5,'color',clrtrace,'linestyle',linestyle);hold on;
        %set(gca,'visible','off');
        xlim([frame.cs_start-16 frame.per_cycle]*fs.ca);ylim(a);box on;
        %saveas(h,[savepathi '\clusnum' num2str(jj,'%02d') '_session' num2str(ii,'%02d') '.eps'],'epsc');
        %close(h);
        t1=mean(B,2);t1=t1(frame.cs_start:frame.us_start-1);t2=b(frame.cs_start:frame.us_start-1);
        %t1=B(frame.cs_start:frame.us_start-1,:);t2=B(frame.cs_start:frame.us_start-1,cIX_iii(ind));
        [ht(ii),pt(ii),~,~]=ttest2(t1,t2,'Vartype','equal');
    end
end
obj = scalebar;
obj.XLen = 2;
obj.YLen = 0.005; %X-Length, 10.
obj.XUnit = 's';            %X-Unit, 'm'.
%obj.Position = [12, -0];  %move the whole SCALE position.
obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
obj.hTextX.FontSize=14;obj.hTextY.FontSize=14;
obj.hTextX.FontSize=16;obj.hTextY.FontSize=16;
obj.hLineX(1).Color=[1 1 1];obj.hLineX(2).Color=[1 1 1];
obj.hLineY(1).Color=[1 1 1];obj.hLineY(2).Color=[1 1 1];
obj.hLineX(1).LineWidth=2;obj.hLineX(2).LineWidth=2;
obj.hLineY(1).LineWidth=2;obj.hLineY(2).LineWidth=2;
obj.hTextX.Color=[1 1 1];obj.hTextY.Color=[1 1 1];
set(gca,'visible','off');

%% acq conditionning change
h_acq=[];
for tt=1:trial.acq_block_num
    x=squeeze(mean(act_all_acq_block_mean(frame.cs_start:frame.us_start-1,tt,:),1));
    y=act_all_acq_block_mean(frame.cs_start-5:frame.cs_start,tt,:);
    h_acq(:,tt)=x>squeeze(mean(y,1)+1*std(y,[],1));
end
num_clust=[2];num_clust=[4,20,10,21,23];%有到增加,
ind=[];zz=1;a=[];%
cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
yli=[-0.02 0.04];kk=1;pt=[];ht=[];ind=[];
for jj=num_clust
    ind=cat(1,ind,find(gIX_iii==jj));
    %correct 2
    %     if jj==4 | jj==7
    %         [ind_c,~,i_ind]=intersect(find(h_acq(:,1)==1),cIX_iii(ind));ind=ind(i_ind);
    % % %     elseif jj==10
    % % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==1),cIX_iii(ind));ind=ind(i_ind);
    % % %     elseif jj==20
    % % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==1),cIX_iii(ind));ind=ind(i_ind);
    % % %     elseif jj==21
    % % %         [ind_c,~,i_ind]=intersect(find(h_acq(:,3)==1),cIX_iii(ind));ind=ind(i_ind);
    %     elseif jj==23
    %         [ind_c,~,i_ind]=intersect(find(h_acq(:,1)==0),cIX_iii(ind));ind=ind(i_ind);
    %         [ind_c,~,i_ind]=intersect(find(h_acq(:,2)==0),cIX_iii(ind));ind=ind(i_ind);
    %         [ind_c,~,i_ind]=intersect(find(h_acq(:,5)==1),cIX_iii(ind));ind=ind(i_ind);
    %      end
    a(zz)=length(ind);  zz=zz+1;
end
rank_ind=find( h_base_adj==0);rank_ind=find(h_base==0 ); rank_ind=find(h_tst==1 );rank_ind=find(h_base==0 & h_tst==1);
rank_ind=find(h_base_adj==0);rank_ind=find(h_base_adj==0 & h_tst_adj==0);rank_ind=find(h_tst_adj==1);rank_ind=find(h_base_adj==1 & h_tst_adj==1);%i_ind=i_ind(500);
[ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(ind));ind=ind(i_ind);
%cur motor related
%figure,hist(motorcorr.sep_CS(cIX_iii(ind),7))
ind(find(motorcorr.sep_CS(cIX_iii(ind),7)>=0.3 ))=[];%& motorcorr.sep_spon(cIX_iii(ind),7)>=0.4
figure,plot(mean(squeeze(act_all_hab_mean(:,:,cIX_iii(ind))),2));hold on;plot(mean(squeeze(act_all_tst_mean(:,:,cIX_iii(ind))),2));hold on;
plot(squeeze(mean(act_all_acq_block_mean(:,:,cIX_iii(ind)),3)));

%ind(find(index_all(cIX_iii(ind))==2 | index_all(cIX_iii(ind))==3 | index_all(cIX_iii(ind))==4 | index_all(cIX_iii(ind))==6))=[];
n=[];
for ii=unique(gIX_iii(ind))'
    n(ii,1)=length(find(gIX_iii(ind)==ii));
end
%% ave. CR
c=7;
figure,
for ii=1:trial.acq_block_num+2
    subplot(1,trial.acq_block_num+2,ii),
    if ii==1
        a=reshape(act_all_hab_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),1,[]);b=squeeze(a(:,1,cIX_iii(ind)))';
    elseif ii==7
        a=reshape(act_all_tst_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),1,[]);b=squeeze(a(:,1,cIX_iii(ind)))';
    else
        a=reshape(act_all_acq_block_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*trial.acq_block_num,[]);a=normalize(a,1);
        a=reshape(a,(frame.us_start-frame.cs_start+c),trial.acq_block_num,[]);b=squeeze(a(:,ii-1,cIX_iii(ind)))';
        %a=act_all_acq_block_mean(frame.cs_start-c:frame.us_start-1,:,:);%
    end
    b=smoothdata(b,1,'gaussian' ,5);
    imagesc(b);hold on;
    line([c c]+1.5,[0 length(ind)],'color','r','linestyle','--','linewidth',2);
    set(gca,'clim',[0.3 1.5]);%set(gca,'visible','off')
end

m=[];sd=[];n=[];
a=reshape(act_all_acq_block_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*trial.acq_block_num,[]);a=normalize(a,1);
a=reshape(a,(frame.us_start-frame.cs_start+c),trial.acq_block_num,[]);
a1=reshape(act_all_hab_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a1=normalize(a1,1);
a1=reshape(a1,(frame.us_start-frame.cs_start+c),1,[]);
a2=reshape(act_all_tst_mean(frame.cs_start-c:frame.us_start-1,:,:),(frame.us_start-frame.cs_start+c)*1,[]);a2=normalize(a2,1);
a2=reshape(a2,(frame.us_start-frame.cs_start+c),1,[]);

for ii=unique(gIX_iii(ind))'
    b=find(gIX_iii(ind)==ii);
    figure,plot(squeeze(a(c+1:15,jj,cIX_iii(ind(b)))))
    for jj=1:5
        m(ii,jj)=mean(sum(a(c+1:15,jj,cIX_iii(ind(b))),1),3);
        sd(ii,jj)=std(sum(a(c+1:15,jj,cIX_iii(ind(b))),1),[],3);
        n(ii)=length(a);
    end
end
figure,plot(squeeze(a(c+1:15,jj,cIX_iii(ind(b)))))
figure,plot(m');figure,plot(mean(m([13,19],:)));hold on;plot(mean(m([11,12,15,17],:)));
m=[];sd=[];n=[];
for ii=unique(gIX_iii(ind))'
    b=find(gIX_iii(ind)==ii); %b=1:length(ind);
    % b=find(gIX_iii(ind)==11 | gIX_iii(ind)==12 | gIX_iii(ind)==15 | gIX_iii(ind)==17);
    %     figure,plot(mean(squeeze(act_all_hab_mean(:,:,cIX_iii(ind(b)))),2));hold on;plot(mean(squeeze(act_all_tst_mean(:,:,cIX_iii(ind(b)))),2));hold on;
    %     plot(squeeze(mean(act_all_acq_block_mean(:,:,cIX_iii(ind(b))),3)));
    n(ii,1)=length(b);
    for jj=1:7
        if jj>=2 && jj<=6
            d=a; m(ii,jj)=mean(mean(d(c+1:15,jj-1,cIX_iii(ind(b))),1),3)-mean(mean(d(3:c,jj-1,cIX_iii(ind(b))),1),3);
            sd(ii,jj)=std(mean(d(c+1:15,jj-1,cIX_iii(ind(b))),1),[],3);
        elseif jj==1
            d=a1;m(ii,jj)=mean(mean(d(c+1:15,1,cIX_iii(ind(b))),1),3)-mean(mean(d(3:c,1,cIX_iii(ind(b))),1),3);
            sd(ii,jj)=std(mean(d(c+1:15,1,cIX_iii(ind(b))),1),[],3);
        else
            d=a2;m(ii,jj)=mean(mean(d(c+1:15,1,cIX_iii(ind(b))),1),3)-mean(mean(d(3:c,1,cIX_iii(ind(b))),1),3);
            sd(ii,jj)=std(mean(d(c+1:15,1,cIX_iii(ind(b))),1),[],3);
        end
    end
end

%% mapback
for kk=[2,3,4,6,12,13];%fish num
    for jj=num_clust;
        a=find(gIX_iii(ind)==jj);a=ind(a);
        %         rank_ind=find(h_base==1 & h_tst==1);
        %         [ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(a));a=a(i_ind);
        
        %a=find(gIX_iii==jj);
        %a=find(gIX_iii==4);a=find(gIX_iii(ind)==10 | gIX_iii(ind)==21);a=find(gIX_iii(ind)==21 | gIX_iii(ind)==23 | gIX_iii(ind)==20);
        
        [i_ind,~,i_a]=intersect(find(index_all==kk),cIX_iii(a)); %i_ind(find(motorcorr.spon(i_ind)>=0.2))=[];figure,hist(motorcorr.sep_CS(i_ind,zz))
        a1=normalize(motorcorr.sep_CS,1);a1=a1(i_ind,7);a2=normalize(motorcorr.sep_spon,1);a2=a2(i_ind,7);
        %c=(a1-a2)./(a2+a1);
        figure,
        %scatter3(X(find(index_all==kk),1),X(find(index_all==kk),2),X(find(index_all==kk),3),8,[1 1 1],'filled');hold on;axis equal;
        scatter3(X(i_ind,1),X(i_ind,2),X(i_ind,3),12,'r','filled');hold on;axis equal;colorbar;colormap('hot');
        
        %         figure,scatter3(X(i_ind,1),X(i_ind,2),X(i_ind,3),12,repmat(cmap(zz,:),length(i_ind),1),'filled');hold on;
        %         zz=zz+1;
        %         axis equal;colorbar;colormap('jet');
        xlim([1 2048]);ylim([1 2048]);grid off;
        set(gca,'visible','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[],'clim',[-1 1]);
        set(gca,'color','k');%axis equal;
        view([0 90]);
    end
end
cmap=[5 3 4 2 1 0]';cmap=[5 2]';
for kk=[2,3,4,6,12,13];%fish num
    figure,   zz=1;
    for jj=num_clust
        a=find(gIX_iii(ind)==jj);a=ind(a);
        %a=find(gIX_iii==jj);
        %         rank_ind=find(h_base==0 & h_tst==1);
        %         [ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(a));a=a(i_ind);
        
        %correct 2
        %[ind_c,~,i_ind]=intersect(find(h_acq(:,zz)==1),cIX_iii(a));a=a(i_ind);
        
        [i_ind,~,i_a]=intersect(find(index_all==kk),cIX_iii(a));
        i_ind(find(motorcorr.sep_CS(i_ind,zz)>=0.3))=[];ind(find(motorcorr.sep_CS(i_ind,7)>=0.3 ))=[];
        %figure,plot_rawtrace_trials(mean(act_all(i_ind,:),1)',[],fs,frame,trial,[],1);
        scatter3(X(i_ind,1),X(i_ind,2),X(i_ind,3),12,repmat(cmap(zz,:),length(i_ind),1),'filled');hold on;
        zz=zz+1;    axis equal;colorbar;colormap('jet');set(gca,'visible','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
        set(gca,'color','k');
        xlim([1 2048]);ylim([1 2048])
        view([-0 90]);
        set(gca,'clim',[-0.5 5]);
    end
end
clrmap_iii=clrmap;clrmap_iii(num_clust,:)=GetColormap('hsv',length(num_clust));
for kk=[2,3,4,6,12,13];
    load(envbatch{kk});
    num_l=find(index_all==kk);
    [~,~,i_ind]=intersect(num_l,cIX_iii(ind));
    for ii=unique(num_clust)
        ii_ind=find(gIX_iii(ind(i_ind))==ii);
        h= DrawTiledPics_zyq_20190530(cIX_iii(ind(i_ind(ii_ind))),gIX_iii(ind(i_ind(ii_ind))),[1:size(A,3)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap_iii);
        saveas(h,[savepath_all '\fish' num2str(kk,'%02d') '_' num2str(ii,'%02d')],'fig');
    end
    %figure,scatter3(X(cIX_iii(ind(i_ind)),1),X(cIX_iii(ind(i_ind)),2),X(cIX_iii(ind(i_ind)),3),'filled');axis equal;
end
%% fraction
%per region
fraction_in_region_in_clust_c=[];ii=1;
for kk=unique(index_all)'
    ind_in_this_fish=(find(index_all==kk));
    cIX_in_this_clus=cIX_iii(ind);%
    [IX,i_ind_in_this_fish,i_cIX]=intersect(ind_in_this_fish' ,cIX_in_this_clus');
    unique(index_all(cIX_in_this_clus(i_cIX)));
    unique(gIX_iii(ind(i_cIX)));
    %e=load(envbatch{kk});
    load([savebatch{kk} '\brain arearegion_mask.mat']);
    %figure,scatter( env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2), env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1))
    %     [~,fraction_in_region_in_clust(:,:,ii),~,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind_in_this_clus(i_cIX)),[1:length(i_cIX)],...
    %         [env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),2),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),1),env_all.supervoxel(cIX_iii(ind_in_this_clus(i_cIX)),3)],clrmap_iii(num_clust,:),true);
    [~,fraction_in_region_in_clust_c(:,:,ii),loc_in_region_cell,~]=get_region_fraction(reg_mask,reg_name,reg_loc,gIX_iii(ind(i_cIX)),cIX_iii(ind(i_cIX)),...
        [env_all.supervoxel(:,2),env_all.supervoxel(:,1),env_all.supervoxel(:,3)],...
        [env_all.supervoxel(ind_in_this_fish,2),env_all.supervoxel(ind_in_this_fish,1),env_all.supervoxel(ind_in_this_fish,3)],clrmap_iii(num_clust,:),false);
    title(['Fish ' num2str(kk,'%02d')]);
    ii=ii+1;
end
l_fraction_in_region_in_clust=mean(fraction_in_region_in_clust_c(:,:,[2,3,4,6]),3,'omitnan');
nl_fraction_in_region_in_clust=mean(fraction_in_region_in_clust_c(:,:,5:end),3,'omitnan');
ii=1;a={};
for jj=num_clust
    a{ii}=[l_fraction_in_region_in_clust(:,jj),nl_fraction_in_region_in_clust(:,jj)];ii=ii+1;
    %plot_fraction_regions(a{ii},Lable,[1 0 0.2;0.2 0 1]);ylim([0 0.1]);
    %title(['Cluster num.' num2str(jj,'%02d')]);
end
x=[];y=[];
for ii=1:size(a,2)
    x=cat(2,x,a{ii}(:,1));y=cat(2,y,a{ii}(:,2));
end
cmap=GetColormap('hsv',length(num_clust));
figure,h=bar(x,'stacked');legend(num2str(num_clust'),'location','bestoutside');
for ii=1:size(h,2)
    h(ii).CData = cmap(ii,:);
end
figure,bar(y,'stacked');legend(num2str(num_clust'),'location','bestoutside');
for ii=1:size(h,2)
    h(ii).CData = cmap(ii,:);
end
for jj=num_clust
    plot_fraction_regions([l_fraction_in_region_in_clust(:,jj),nl_fraction_in_region_in_clust(:,jj)],Lable,[0.5 0.5 0.5;0.2 0 1]);ylim([0 0.2]);
    %plot_fraction_regions([l_fraction_in_region_in_clust(:,jj)-nl_fraction_in_region_in_clust(:,jj)]./[l_fraction_in_region_in_clust(:,jj)+nl_fraction_in_region_in_clust(:,jj)],Lable,[0.5 0.5 0.5;0.2 0 1]);ylim([-0.5 1.5]);
    title(['Cluster num.' num2str(jj,'%02d')]);
end

%% find spacial cluster
num_clust=[4,10,20,21,23];ind=[];zz=1;
for jj=num_clust
    ind=cat(1,ind,find(gIX_iii==jj));
    %     % correct 1
    %     rank_ind=find(h_base==0 & h_tst==1);%i_ind=i_ind(500);
    %     [ind_c,~,i_ind]=intersect(rank_ind,cIX_iii(ind));ind=ind(i_ind);
    %     %correct 2
    %     [ind_c,~,i_ind]=intersect(find(h_acq(:,zz)==1),cIX_iii(ind));ind=ind(i_ind);zz=zz+1;
end

figure,plot_rawtrace_trials(mean(act_all(cIX_iii(ind),:),1)',[],fs,frame,trial,[],1);

c=[];
for ii=1:2
    if ii==1
        B=squeeze(act_all_hab_mean);B=normalize(B,1);
    else
        B=squeeze(act_all_tst_mean);B=normalize(B,1);
    end
    b=mean(B(frame.cs_start:frame.cs_end,ind),2);
end
i_ind=find(mean(c,2)>0.6);%i_ind=i_ind(1);%figure,imagesc(corr(act_all(i_ind,:)))
figure,plot_rawtrace_trials(mean(act_all(i_ind,:),1)',[],fs,frame,trial,[],1);
