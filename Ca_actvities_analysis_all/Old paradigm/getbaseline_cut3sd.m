%�ж�CSǰ��û���Է���Ӧ�Լ��쳣��Ⱥֵ,����У�ȥ��
%aΪ�������������еľ���indcutΪҪ����ķ�Χ����������
%�ж���������Ϊ��ʼ��>1sd��thr2������ߵ�>2sd��thr1��,�Ƿ�Ӧ�Ļ�������
function base=getbaseline_cut3sd(aa,indcut)
%a=activities_total{ii}(:,jj);indcut=1:frame.cs_start-1;
%aa=reshape(a,frame.per_cycle,trial.total)';

bb=aa(:,indcut);b=smoothdata(bb,2,'movmean',3);
b=sort(bb,2);b=b(:,1:fix(length(indcut)*0.95));%ȡ������С��һ����Ϊ�жϱ�׼
m=mean(b,2);sd=std(b,[],2);
thres1=m+3*sd;thres1=kron(thres1,ones(1,length(indcut)));
thres2=m+2*sd;thres2=kron(thres2,ones(1,length(indcut)));

% figure,plot(reshape(aa',1,[])','b');hold on;plot(reshape(kron(m+3*sd,ones(1,50))',1,[])');hold on;
% plot(reshape(smoothdata(aa,2,'movmean',3)',1,[])')

b=reshape(bb',1,[]);thr1=reshape(thres1',1,[]);thr2=reshape(thres2',1,[]);
m=reshape(kron(m,ones(1,length(indcut)))',1,[]);
ind=find(abs(b)>=thr1);
b(ind)=m(ind);%����thr1��ȫ����Ҫ
%ȥǰ3����ͺ������㣬��ֹ����
if ~isempty(find(ind<=3))
    %b(ind(find(ind<=3)))=m(ind(find(ind<=3)));
    ind(find(ind<=3))=[];
end
if ~isempty(find(ind>=length(b)-2))
    %b(ind(find(ind>=length(b)-2)))=m(ind(find(ind>=length(b)-2)));
    ind(find(ind>=length(b)-2))=[];
end
%��Ϊǰ������е�>thr2�������Ƿ�Ӧ������ȥ�������ֻ�ǵ����㣬����ֻ���������Ƕ��ˣ�
if ~isempty(ind)
    b(ind(find(b(ind-2)>thr2(ind-2)))-2)=m(ind(find(b(ind-2)>thr2(ind-2)))-2);
    b(ind(find(b(ind-1)>thr2(ind-1)))-1)=m(ind(find(b(ind-1)>thr2(ind-1)))-1);
    b(ind(find(b(ind+1)>thr2(ind+1)))+1)=m(ind(find(b(ind+1)>thr2(ind+1)))+1);
    b(ind(find(b(ind+2)>thr2(ind+2)))+2)=m(ind(find(b(ind+2)>thr2(ind+2)))+2);
end
bb=reshape(b',length(indcut),size(aa,1))';
aa1=aa;aa1(:,indcut)=bb;aaa=reshape(aa1',1,[])';

base.output_cut3sd=aa1;
base.input_ind_cut3d=bb;
base.m=mean(bb,2);base.sd=std(bb,[],2);
base.rep_m=reshape(kron(base.m,ones(1,size(aa,2)))',1,[])';
base.rep_sd=reshape(kron(base.sd,ones(1,size(aa,2)))',1,[])';

% figure,plot(reshape(aa',1,[])','b');hold on;plot(aaa,'r');hold on;
% baseline1=mean(aa(:,indcut),2)+1*std(aa(:,indcut),[],2);
% baseline1=kron(baseline1,ones(1,size(aa,2)));
% baseline1=reshape(baseline1',1,[]);
% plot(baseline1,'c');hold on;
% baseline2=base.m+1*base.sd;
% baseline2=kron(baseline2,ones(1,size(aa,2)));
% baseline2=reshape(baseline2',1,[]);
% plot(baseline2,'m');hold on;




