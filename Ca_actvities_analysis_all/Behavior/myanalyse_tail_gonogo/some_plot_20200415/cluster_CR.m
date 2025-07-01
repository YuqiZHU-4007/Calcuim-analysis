%% load movement parament (theta & duration)
%load
path=uigetdir('H:\1.Test US\2.Tail free¡ª¡ªData from 117');
seg={'20200618','20200624','20200625','20200626','20200627'};
n=[4,5,4,6,5];kk=1;
Path=string;Name=string;
for ii=1:length(n)
    for jj=1:n(ii)
        name=fullfile(path,seg{ii},['fish' num2str(jj)],'Results_of_alltheta.mat');
        Path(kk,1)=name;
        Name(kk,1)=fullfile(seg{ii},['fish' num2str(jj)]);
        kk=kk+1;
    end  
end
load(Path(2));name=fieldnames(Results_of_alltheta);
Data=struct;P=[];alltheta_all=nan(500,36,length(Path));
for ii=1:length(name)
     Data=setfield(Data,name{ii},[]);
end
for zz=1:length(Path)
    if exist(Path(zz))
        load(Path(zz));ind=find(strcmp(Path,Path(zz)));
        for ii=1:length(name)
            if min(size(getfield(Results_of_alltheta,name{ii})))<=1
                b=getfield(Data,name{ii})';a=getfield(Results_of_alltheta,name{ii});if ~isempty(a) b(ind,1:length(a))=a;end
            elseif min(size(getfield(Results_of_alltheta,name{ii})))>1
                b=getfield(Data,name{ii})';
                for jj=1:trial.acq_block_trial
                    bb=getfield(Results_of_alltheta,name{ii});
                    b{ind,jj}=bb(:,jj,:);
                end
            elseif min(size(getfield(Results_of_alltheta,name{ii})))>1 && strcmp(name{ii},'spon_env_loc')
                b=getfield(Data,name{ii})';
                for jj=1:size(getfield(Results_of_alltheta,name{ii}),1)
                    bb=getfield(Results_of_alltheta,name{ii});
                    b{ind,jj}=bb(jj,:,:);
                end
            end
            Data=setfield(Data,name{ii},b');
        end
        for hh=3
            switch hh
                case 1
                    ind_cs=ind_CS_hab;
                case 2
                    ind_cs=ind_CS_tst;
                case 3
                    ind_cs=ind_CS;
            end
            aa=(repmat(T_non_spon(env.env_loc(:,1)),length(ind_cs(:,1)),1)-repmat(ind_cs(:,1),1,length(env.env_loc(:,1))));
            ind=find(aa>=-0.1*60 & abs(aa)<=(frameb.cs_end-frameb.cs_start)+1-0.2*60);[~,ind]=ind2sub([size(aa,1),size(aa,2)],ind);
            env_loc_cs=[env.env_loc(ind,:),env.env_end(ind,:)];env_num=length(ind);
            [p,x]=calcula_bout_para(alltheta_n,T_non_spon(env_loc_cs(:,1)),T_non_spon(env_loc_cs(:,2)));
            [C,ia,~]=unique(fix(env_loc_cs(:,1)/frameb.per_cycle)+1);
            for ii=C'
                ind=find(fix(env_loc_cs(:,1)/frameb.per_cycle)+1==ii);
                %P(1:length(fieldnames(p)),ii,zz)=mean(squeeze(cell2mat(squeeze(struct2cell(p(ind))))),2);
                P(1,ii,zz)=max([p(ind).max_amp]);
                P(2,ii,zz)=mean([p(ind).mean_amp]);
                P(3,ii,zz)=max([p(ind).max_tail_change]);
                P(4,ii,zz)=nanmean([p(ind).mean_tail_change]);
                P(5,ii,zz)=sum([p(ind).sum_tail_change]);
                P(6,ii,zz)=sum([p(ind).dur]);
                P(7,ii,zz)=sum([p(ind).beat_num]);
                P(8,ii,zz)=sum([p(ind).beat_num_2]);
                P(length(fieldnames(p))+1,ii,zz)=length(ind);
                a=T_non_spon(env_loc_cs(ind(1),1)):T_non_spon(env_loc_cs(ind(end),2));
                alltheta_all(1:length(a),ii,zz)=alltheta_n(a);
            end
        end
    end
end
%% PCA
size(P)
%reshape
Para=reshape(P,size(P,1),[],size(P,4));Para=reshape(Para,size(P,1),[],1);
Para=Para';
alltheta_all_rep=reshape(alltheta_all,size(alltheta_all,1),[],size(alltheta_all,4));alltheta_all_rep=reshape(alltheta_all,size(alltheta_all,1),[],1);
alltheta_all_rep(:,find(Para(:,9)==0))=[];Para(find(Para(:,9)==0),:)=[];alltheta_all_rep=alltheta_all_rep';

%Para=[Para;A';B']';
Para=normalize(Para,'zscore');
[coeff,score,latent,~,explained] = pca(Para);
explained_var = cumsum(latent)/sum(latent)
figure('position',[100,100,1200,430]),subplot(1,2,1),plot(explained_var,'-k', 'Marker','o','MarkerSize',5,'linewidth',2);hold on;
bar(explained/100,'FaceColor',[0.5 0.5 0.5]);
plot(min(find(explained_var>0.9)),explained_var(min(find(explained_var>0.9))),'-r', 'Marker','o','MarkerSize',5,'linewidth',2);hold off;
xlim([0 25]);xlabel('# PC','fontsize',15);ylabel('variance explained (%)','fontsize',15);
subplot(1,2,2),imagesc(coeff);colorbar;xlabel('# PC','fontsize',15);ylabel('Parameters','fontsize',15);

figure,scatter3(score(:,1),score(:,2),score(:,3));hold on; 
xlabel('PC1');ylabel('PC2');zlabel('PC3');
figure,
kk=1;
for ii=1:4
    for jj=1:4
        subplot(4,4,kk),
        scatter(score(:,ii),score(:,jj));hold on;
        xlabel(['PC' num2str(ii)]);ylabel(['PC' num2str(jj)]);
        kk=kk+1;
    end
end
%% clustering
index_UR_2=1:size(Para,1);
PC_num=2:4;
X=score(:,PC_num);
range=2:10;Silh_mean=[];
for i = 1:length(range),
    disp(['k-means k=' num2str(range(i))]);
    [ix,~,~]= kmeans(X,range(i),'Distance','sqeuclidean','Replicates',10,'Start','plus');
    silh = silhouette(X,ix,'sqeuclidean');
    Silh_mean(i) = mean(silh);
    subplot(4,4,i),col=hsv(range(i));
    scatter3(X(:,1),X(:,2),X(:,3),12,col(ix,:),'filled');
    xlabel(['PC' num2str(PC_num(1))]);ylabel(['PC' num2str(PC_num(2))]);zlabel(['PC' num2str(PC_num(3))]);
end
figure;plot(range,Silh_mean,'o-','color',[0.5 0.5 0.5]);

[~,ix] = max(Silh_mean);n=9;
col=hsv(n);
[idx_kmeans,C,sumd,D] = kmeans(X,n,'Distance','sqeuclidean','Replicates',10,'Start','plus');
figure,
plot3(C(:,1),C(:,2),C(:,3),'kx','Markersize',15);hold on;
scatter3(X(:,1),X(:,2),X(:,3),12,col(idx_kmeans,:),'filled');
xlabel(['PC' num2str(PC_num(1))]);ylabel(['PC' num2str(PC_num(2))]);zlabel(['PC' num2str(PC_num(3))]);
figure,
kk=1;
for ii=1:size(X,2)
    for jj=1:size(X,2)
        subplot(size(X,2),size(X,2),kk),
        gscatter(X(:,ii),X(:,jj),idx_kmeans);hold on;
        xlabel(['PC' num2str(PC_num(ii))]);ylabel(['PC' num2str(PC_num(jj))]);
        kk=kk+1;legend off;
    end
end
legend('Cluster 1','Cluster 2','Cluster 3','Labeled')
figure,
for ii=1:n
    a=1:size(alltheta_all_rep,2);
    b= mean(abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:)),'omitnan');
    if length(find(idx_kmeans==ii))==1
       plot(a,b,'color',col(ii,:));hold on;
    else
    errBar=std(abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:)),'omitnan')/sqrt(length(find(idx_kmeans==ii))-1);
    shadedErrorBar(a,b,errBar,'lineprops',col(ii,:));hold on;
    end
    
end
%mapback
Para=reshape(P,size(P,1),[],size(P,4));Para=reshape(Para,size(P,1),[],1);
ind_clus=find(Para(:,9)~=0);
for ii=1
    index_UR_2(find(idx_kmeans==ii))
end
%% svm