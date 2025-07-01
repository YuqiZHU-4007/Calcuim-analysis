%input
nCells_total ; 
prct_const;
reg_thres;
M;
%% regressor
behavior=zeros(5,length(time_be));behavior(1,:)=abs(y);
[regressors,~,~,~] = GetMotorRegressor(behavior);
output=[regressors(1,1).im];
[corr,~] = MotorSourceCorrelation(M,output,[]);
[~,IX] = sort(corr,'descend');
topN = round(prct_const/100 * nCells_total);
x0 = corr(IX(topN));
ind_output = (find(corr>max(reg_thres,x0)))';
output=mean(M(ind_output,:),1);
