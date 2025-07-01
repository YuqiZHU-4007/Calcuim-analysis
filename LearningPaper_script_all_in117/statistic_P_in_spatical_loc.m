function [P_z,H_z,P_t,H_t,Dr_ij,Dn_ij,Dr_ij_projection,Dn_ij_projection]=statistic_P_in_spatical_loc(loc,all_loc)
%loc:1*n cell,每个cell中为需要计算p的神经元坐标xyz，不同cell为不同个体
%all_loc：对应每个个体所有神经元坐标
Dr_ij={}; Dn_ij={};Dn_ij_parfor={};N=10000;
for fish_i=1:length(loc)
    all_n_in_fish_i=loc{fish_i};
    for neuron_in_fish_i=1:length(all_n_in_fish_i)
        kk=1;
        for fish_j=setdiff(1:length(loc),fish_i)
            other_in_fish_j=loc{fish_j};
            Dr_ij{fish_i}(neuron_in_fish_i,kk)=min(pdist2(all_n_in_fish_i(neuron_in_fish_i,:),other_in_fish_j));kk=kk+1;
        end
    end
end

parfor nsample=1:N
    for fish_i=1:length(loc)
        all_n_in_fish_i=loc{fish_i};
        for neuron_in_fish_i=1:length(all_n_in_fish_i)
            kk=1;
            for fish_j=setdiff(1:length(loc),fish_i)
                rand_id=randi(length(all_loc{fish_j}),length(loc{fish_j}),1);
                random_in_fish_j=all_loc{fish_j}(rand_id,:);
                Dn_ij_parfor{nsample}{fish_i}(neuron_in_fish_i,kk)=min(pdist2(all_n_in_fish_i(neuron_in_fish_i,:),random_in_fish_j));kk=kk+1;
            end
        end
    end
end
for fish_i=1:length(loc)
    aa=[];
    for ii=1:N
        aa(:,:,ii)=Dn_ij_parfor{ii}{fish_i};
    end
    Dn_ij{fish_i}=aa;
end

%% projection fucntion of D
percentage=6/8;A=0;
for fish_i=1:length(loc)
    all_n_in_fish_i=loc{fish_i};
    A=max(A,length(all_n_in_fish_i));
end
Dr_ij_projection=nan(length(loc),A);Dn_ij_projection=nan(length(loc),A,N);
for fish_i=1:length(loc)
    a=sort(Dr_ij{fish_i},2,'ascend');
    Dr_ij_projection(fish_i,1:size(a,1))=mean(a(:,1:ceil(size(a,2)*percentage)),2); 
    a=sort(Dn_ij{fish_i},2,'ascend');
    Dn_ij_projection(fish_i,1:size(a,1),:)=squeeze(mean(a(:,1:ceil(size(a,2)*percentage),:),2));
end
%% 重新统计P：ztest 每个neuron的null distribution来自于这个neruon的n sampling
alpha=0.05;
P_z=nan(length(loc),A);H_z=nan(length(loc),A);%单尾，h1：u<u0(Dr<Dn)
P_t=nan(length(loc),A);H_t=nan(length(loc),A);
for fish_i=1:length(loc)
    all_n_in_fish_i=loc{fish_i};
    for neuron_in_fish_i=1:length(all_n_in_fish_i)
        a=squeeze(Dn_ij_projection(fish_i,neuron_in_fish_i,:));m=mean(a,'omitnan');s=std(a,'omitnan');n=length( Dr_ij{fish_i}(neuron_in_fish_i,:));
        t=(Dr_ij_projection(fish_i,neuron_in_fish_i)-m)./(s/sqrt(n));
        p=normcdf(t);       
        t_a=norminv(alpha);if t<-abs(t_a);  h=true;else h=false;end
        P_z(fish_i,neuron_in_fish_i)=p;
        H_z(fish_i,neuron_in_fish_i)=h;
    end
end
a=find(P_z<alpha);
b=find(H_z==1);

if false
    %% 统计
    alpha=0.05;
    P_z=nan(length(loc),A);H_z=nan(length(loc),A);%单尾，h1：u<u0(Dr<Dn)
    for fish_i=1:length(loc)
        all_n_in_fish_i=loc{fish_i};
        for neuron_in_fish_i=1:length(all_n_in_fish_i)
            a=squeeze(Dn_ij{fish_i}(neuron_in_fish_i,:,:));a=reshape(a,1,[]);
            m=mean(a);s=std(a);n=length( Dr_ij{fish_i}(neuron_in_fish_i,:));
            t=(mean( Dr_ij{fish_i}(neuron_in_fish_i,:),2)-m)./(s/sqrt(n));p=normcdf(t);
            t_a=norminv(alpha);if t<-abs(t_a);  h=true;else h=false;end
            P_z(fish_i,neuron_in_fish_i)=p;
            H_z(fish_i,neuron_in_fish_i)=h;
        end
    end
    a=find(P_z<alpha);
    b=find(H_z==1);
    
    P_t=nan(length(loc),A);H_t=nan(length(loc),A);%单尾，h1：u<u0(Dr<Dn)
    for fish_i=1:length(loc)
        all_n_in_fish_i=loc{fish_i};
        for neuron_in_fish_i=1:length(all_n_in_fish_i)
            a=squeeze(Dn_ij{fish_i}(neuron_in_fish_i,:,:));a=reshape(a,1,[]);
            m=mean(a);s=std(a);n=length( Dr_ij{fish_i}(neuron_in_fish_i,:));
            t=(mean( Dr_ij{fish_i}(neuron_in_fish_i,:),2)-m)./(s/sqrt(n));
            t_a=tinv(alpha,n-1);
            if t<-abs(t_a); h=true;else h=false;end
            p=tcdf(t,n-1);
            P_t(fish_i,neuron_in_fish_i)=p;
            H_t(fish_i,neuron_in_fish_i)=h;
        end
    end
    % a=find(P_z<alpha);
    % b=find(H_z==1);
end
% for nsample=1:N
%     for fish_i=1:length(loc)
%         all_n_in_fish_i=loc{fish_i};
%         for n_in_fish_i=1:length(all_n_in_fish_i)
%             kk=1;
%             for fish_j=setdiff(1:length(loc),fish_i)
%                 rand_id=randi(length(all_loc{fish_j}),length(loc{fish_j}),1);
%                 random_in_fish_j=all_loc{fish_j}(rand_id,:);
%                  Dn_ij{fish_i}(n_in_fish_i,kk,nsample)=min(pdist2(all_n_in_fish_i(n_in_fish_i,:),random_in_fish_j));kk=kk+1;
%             end
%         end
%     end
% end