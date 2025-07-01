function [P_z,H_z,Dr_ij_projection,Dn_ij_projection]=adjust_P_spatial_loc(loc,Dr_ij,Dn_ij)
N=size(Dn_ij{1},3);
%% projection fucntion of D
percentage=0.6;A=0;
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