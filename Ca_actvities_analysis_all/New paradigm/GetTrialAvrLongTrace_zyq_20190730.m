function [C_trialAvr,C_trialRes,C_score,C_d2var_perstim,C_p] = GetTrialAvrLongTrace_zyq_20190730(C,periods)

period = periods;
C_3D_0 = reshape(C,size(C,1),period,[]);
%     C_3D = zscore(C_3D_0,0,2);
%     C_d2var_perstim = nanmean(nanstd(C_3D,0,3),2);
%     C_score = C_d2var_perstim;

C_period = mean(C_3D_0,3);%prctile(C_3D_0,20,3);%mean(C_3D_0,3);
nPeriods = round(size(C,2)/period);
C_trialAvr = repmat(C_period,1,nPeriods);

C_trialRes = C-C_trialAvr;
C_3D_tRes = reshape(C_trialRes,size(C,1),period,[]);
%     C_3D_tRes = zscore(C_3D_tRes_0,0,2);
d2var = (nanstd(C_3D_tRes,0,3)).^2;
C_score = nanmean(d2var,2);

C_p = C;
C_d2var_perstim = [];