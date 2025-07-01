%
[frame,trial]=setpara_testUS(3);
[stimUS, ind_US,T_spon,T_US,l_spon,l_us]=getUSstim(trial,frame);%figure,plot( stimUS);

n=char(Path);Name=string(n(:,end-14:end-1))';
unresponsive=[63 70 77 78];N=Name(unresponsive);
type1=[56:59 62 68 74 76 79 80 83  78 ];N=Name(type1);%4 week
type2= [65:67 71 73];N=Name(type2);%3 week
b=Data.M_num_US(1:8,type2);b=Data.M_num_US(1:8,unresponsive);

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
type=type2;
b=Data.avg_env_dur_us(1:8,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),8,1))./(b+repmat(b(1,:),8,1));

b=Data.Theta_avg_US_n(1:8,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),8,1))./(b+repmat(b(1,:),8,1));

b=Data.M_num_spon(1:9,type);
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

b=Data.Theta_avg_spon_n(1:9,type);b(isnan(b))=0;
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

b=Data.Theta_US(:,type);b=cell2mat(b);b=reshape(b,trial.acq_block_number,trial.acq_block_trial,[]);
figure,hist3(b)
xlabel('MPG')
ylabel('Weight')

x=[0	0	0	0.097555309152363	0	0	0.952343814981903	0	0.0363442298556288	0.203164126191337];
y=[0	0	0	0.196365300329236	0.221049474968585	0.0635301393834551	0	0.0017394766026371	0	0];
[h,p]=ttest2(x,y,'tail','right')