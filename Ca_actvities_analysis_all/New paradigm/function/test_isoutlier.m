function test_isoutlier(x,fs,frame,trial)
xx=reshape(x,1,[]);

TF=isoutlier(x);
ind = find(TF);
x1 = filloutliers(x,'next');
x1=reshape(x1,1,[]);

figure,plot_rawtrace_trials([x1;xx]',[],fs,frame,trial,[],1);

figure,plot(xx,'b');hold on
plot(x1,'--r');hold on
legend('Noisy Data with Outlier','Noisy Data with Filled Outlier')