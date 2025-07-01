function [start_index,y,m,sd]=findonset(data,basewin,win)
% 就是利用峰值点之前差分值大于0，峰值点之后差分小于0，来判断峰值点位置。
y=smoothdata(data,1);
start_index=nan(1,size(data,2));m=nan(1,size(data,2));sd=nan(1,size(data,2));
%figure,plot(data,'r');hold on;plot(y,'b');xlim([5,70])
m=mean(y(basewin,:),1,'omitnan');sd=std(y(basewin,:),[],1,'omitnan');
c=3;
for i=1:size(y,2)
    b=y(:,i);
    a=find(b>m(i)+1*sd(i) & [b(2:end);b(end)]>m(i)+2*sd(i)  & [b(3:end);b(end);b(end)]>m(i)+2*sd(i));
    if ~isempty(find(a>min(win),1))
    start_index(1,i)= a(find(a>min(win),1));
    end
end
% figure,plot(y(:,4),'b--');hold on;plot(data(:,5),'r');xlim([5,70]);
% scatter(start_index(4),y(start_index(4),4))
end