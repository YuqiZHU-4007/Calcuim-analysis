%先去高值点求阈值，再判断反应，把CS之前的反应去掉
function compute_dfdf_preCS
a=activities_total{ii};
aa=reshape(a(:,jj),frame.per_cycle,trial.total)';

%[~,b]=compute_dfdf(aa,11,0.5);bb=b(:,1:frame.cs_start-1);
b=aa(:,1:frame.cs_start-1);bb=sort(b,2);bb=bb(:,1:fix(frame.cs_start*0.90));

%figure,plot(aa(2,:));hold on;plot(bb(2,:));hold on;
%line([1 50],[mean(aa(2,:)) mean(aa(2,:))]);hold on; line([1 50],[mean(aa(2,1:25)) mean(aa(2,1:25))],'color','r');

b(:,frame.cs_start:frame.per_cycle )= repmat(mean(bb,2),1,length(frame.cs_start:frame.per_cycle))+3*repmat(std(bb,[],2),1,length(frame.cs_start:frame.per_cycle));
m=repmat(mean(bb,2),1,size(aa,2));sd=repmat(std(bb,[],2),1,size(aa,2));mm=reshape(m',size(aa,1)*size(aa,2),1);
ind=zeros(size(aa));ind=abs(b-m)-2*sd>0; ind(:,frame.cs_start:frame.per_cycle)=0;ind=reshape(ind',size(aa,1)*size(aa,2),1);%bb=reshape(b',size(b,1)*size(b,2),1);
if ~isempty(find(ind==1))
    a(find(ind==1),jj)=mm(find(ind==1));%interp1(1:length(a(:,jj)),a(:,jj),find(ind==1),'linear');
   figure,plot(activities_total{ii}(:,jj));hold on;plot(a(:,jj),'r');
end
aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
activities_preCS_baseline{ii,1}(:,jj)=reshape(repmat(mean(aa(:,1:frame.cs_start-1),2),1,size(aa,2))',frame.per_cycle*trial.total,1);
activities_preCS_dfdf{ii,1}(:,jj)=(activities_total{ii}(:,jj)-activities_preCS_baseline{ii}(:,jj))./activities_preCS_baseline{ii}(:,jj);
figure,plot(activities_total{ii}(:,jj));hold on;plot(activities_preCS_baseline{ii}(:,jj),'r');