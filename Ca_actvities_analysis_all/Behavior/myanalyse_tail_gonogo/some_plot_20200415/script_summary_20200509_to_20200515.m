%
n=char(Path);Name=string(n(:,end-14:end-1))';

type1=[32 38 41 43 46 49 51 52 53];
type2= [31 34 35 39 44 45 50];
type3=[36 37 47 48];
unchanged=[38 41 46 52];Name(unchanged)
unsure=[32 43 49 51];Name(unsure)
changed=[53 31 34 35 39 44 45 50];Name(changed)
type_all=[31 32 34:39 41 43:53];[type1,type2,type3];
b=Data.M_num_US(1:8,[type1,type2,type3]);
b(isnan(b))=0;
figure,
for ii=1:trial.acq_block_number
    bbb=[];
    for jj=0:trial.acq_block_trial
        bbb(jj+1)=length(find(b(ii,:)==jj))/length(b(ii,:));
    end
    subplot(1,trial.acq_block_number,ii),pie(bbb);
    title(['Block.' num2str(ii)],'fontsize',14);
end
labels = {'0','1','2','3'};
legend(labels)

type=unchanged;
type=changed;
type=[type1,type2];
type=unsure;
type=type1;
b=Data.avg_env_dur_us(1:8,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),8,1))./(b+repmat(b(1,:),8,1));

b=Data.Theta_avg_US_n(1:8,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),8,1))./(b+repmat(b(1,:),8,1));

b=Data.M_num_spon(1:9,type);
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

b=Data.Theta_avg_spon_n(1:9,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

x=[6.6	2.2	2.4	1	5.2	0.4	3.2	4	4	4	8.2	0.2	2.4	5.2	5.4	3];
y=[4.2	3.4	0.6	1.4	5.4	3.4	1.2	1.2	4.6	0	4.2	0.4	0	4.4	0.2	0];
[h,p]=ttest2(x,y,'tail','right')