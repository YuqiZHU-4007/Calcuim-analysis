function [trail_mov_number,trail_nomov_number,trail_false_number,sti_tdur,base_tdur]=find_mov_trailnum(t_dur_cell,sti_start,duration,trailnum)
k1=1;k2=1;k3=1;
trail_mov_number=zeros;trail_nomov_number=zeros;trail_false_number=zeros;
for i=1:trailnum
    a=t_dur_cell{i,:};%%%cell
    a(find(a(:,1)==0),:)=[];
    ll=find(a(:,1)>=sti_start & a(:,1)<=(sti_start+duration));
    sti_tdur(1:length(ll),:,i)=a(ll,:);
    ll2=find(a(:,1)<sti_start);
    base_tdur(1:length(ll2),:,i)=a(ll2,:);
    ll3=find(a(:,1)>(sti_start+duration)); 
    if length(ll)~=0
        trail_mov_number(k1)=i;
        k1=k1+1;
    elseif length(ll)==0 && (length(ll2)~=0 || length(ll3)~=0)
        trail_nomov_number(k2)=i;
        k2=k2+1;
    elseif length(ll)==0 && length(ll2)==0 && length(ll3)==0
        trail_false_number(k3)=i;k3=k3+1;
    end
end