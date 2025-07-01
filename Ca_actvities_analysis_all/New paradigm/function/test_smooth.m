function [x1,x2,x3,error]=test_smooth(x,fs,frame,trial,polyfit_win,ind)
xx=reshape(x,1,[]);
x1=smoothdata(x,'sgolay');
x1=reshape(x1,1,[]);

x2=smoothdata(x,'rlowess',3);
x2=reshape(x2,1,[]);

x3=smoothdata(x,'movmean',3);
x3=reshape(x3,1,[]);


% figure,plot_rawtrace_trials([x1;xx]',[],fs,frame,trial,[],1);
% figure,plot_rawtrace_trials([x2;xx]',[],fs,frame,trial,[],1);
% figure,plot_rawtrace_trials([x3;xx]',[],fs,frame,trial,[],1);

% figure,plot(xx,'b');hold on
% plot(x1,'--g');hold on
% plot(x3,'--r');hold on
% plot(x2,'--k');hold on
% legend('Data','sgolay','rlowess_3','movmean_3');
figure_alble=true;
if figure_alble
figure,
end

X=[xx;x1;x2;x3];
for ii=1:4
    a=X(ii,:);
    [p,S] = polyfit(polyfit_win,a(polyfit_win),3);
    [lower_integrate,delta]=polyval(p,[polyfit_win ind],S);%
    %lower_integrate=p(1)*[polyfit_win ind]+p(2);
    %error(:,ii)=sqrt(sum((xx(ind)-lower_integrate(length(polyfit_win)+1:end)).^2));
    error(:,ii)=mean((xx(ind)-lower_integrate(length(polyfit_win)+1:end)));
    
    if figure_alble
    subplot(2,2,ii),
    plot(1:length(xx),xx,'b','linewidth',2);hold on
    plot(1:length(xx),a,'--k','linewidth',2);hold on
    hold on
    plot([polyfit_win ind],lower_integrate,'r-','linewidth',2)
    plot([polyfit_win ind],lower_integrate+2*delta,'m--',...
        [polyfit_win ind],lower_integrate-2*delta,'m--','linewidth',2)
    title('Linear Fit of Data with 95% Prediction Interval')
    legend('Data-smooth','Data','Linear Fit','95% Prediction Interval');
    text(20,0.005,num2str(error(:,ii)),'fontsize',10);
    ylim([min(xx)-0.5*(max(xx)-min(xx)) max(xx)+0.5*(max(xx)-min(xx))])
    ylim([-0.03 0.06])
    end
end

