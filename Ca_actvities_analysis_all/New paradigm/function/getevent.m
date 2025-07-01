function event=getevent(input,trial,frame)

aa=reshape(input,frame.per_cycle,trial.total)';
base=getbaseline_cut3sd_strict_for_new_paradigm(aa,1:frame.cs_start-1,trial);
m=base.rep_m;sd=base.rep_sd;
% m=mean(aa(:,1:frame.cs_start-1),2);m=kron(m,ones(frame.per_cycle,1));
% sd=std(aa(:,1:frame.cs_start-1),[],2);sd=kron(sd,ones(frame.per_cycle,1));
baseline=m+2*sd;maxthres=m+3*sd;
rasing_width=[1 100];
decling_width=[4 100];%decling_width+rasing_width-1<8
event=getevent_Trianglefitting_strict_for_new_paradigm(input,baseline,maxthres,rasing_width,decling_width,frame,trial,base);