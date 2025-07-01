close all;
set(0,'defaultfigurecolor','w');
shake_cs_spon123mean_std={};shake_cs_mean_std={};
%1:learner;
%2:nonlearner
for ii=1:4
    switch ii
        case 1
            x1=shake_spon_flash_dimming_learner_all;
            x2=shake_cs_flash_dimming_learner_all;
            shake_cs_mean_std{ii} = shake_individual_bin(x2,trial);
        case 2
            x1=shake_spon_flash_dimming_nonlearner_all;
            x2=shake_cs_flash_dimming_nonlearner_all;
            shake_cs_mean_std{ii} = shake_individual_bin(x2,trial);
        case 3
            x1=shake_spon_flash_dimming_learner_all;
            x2=[];
        case 4
            x1=shake_spon_flash_dimming_nonlearner_all;
            x2=[];
    end
    xx=cat(3,x1,x2);
    xl=min(sum(xx,3),1,'includenan')/size(xx,3);
    shake_cs_spon123mean_std{ii} = shake_individual_bin(xl,trial);
end

a1=shake_cs_spon123mean_std{1};
a2=shake_cs_spon123mean_std{2};
yy1=shake_cs_mean_std{1};
yy2=shake_cs_mean_std{2};
x1=a1(1,:,1);x2=a2(1,:,1);xx=[x1,x2];
y1=yy1(7,:,1);y2=yy2(7,:,1);yy=[y1,y2];

p1 = polyfit(x1,y1,1);
yfit1 = polyval(p1,x1);
yresid = y1 - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(y1)-1) * var(y1);
rsq1 = 1 - SSresid/SStotal

p2 = polyfit(x2,y2,1);
yfit2 = polyval(p2,x2);
yresid = y2 - yfit2;
SSresid = sum(yresid.^2);
SStotal = (length(y2)-1) * var(y2);
rsq2 = 1 - SSresid/SStotal

p = polyfit(xx,yy,1);
yfit = polyval(p,xx);
yresid = yy - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(yy)-1) * var(yy);
rsq = 1 - SSresid/SStotal

figure,
s1=scatter(x1,y1,'r','filled');hold on;
s2=scatter(x2,y2,'b','filled');hold on;
l1=plot(x1,yfit1,'Color','r','linewidth',2);hold on;
l2=plot(x2,yfit2,'Color','b','linewidth',2);hold on;
legend([s1,s2,l1,l2],{'learner','non-learner',...
    ['linear fit: y=' num2str(p1(1)) '*x+' num2str(p1(2))],['linear fit: y=' num2str(p2(1)) '*x+' num2str(p2(2))]});
text(0.75,0.8,['R^2-learner=' num2str(rsq1)]);
text(0.75,0.70,['R^2-nonlearner=' num2str(rsq2)]);
axis equal;xlim([0 1]);ylim([0 1]);
xlabel('Hab. shaking performance');
ylabel('Tst. learning performance');
set(gca,'fontsize',12);

figure,
s1=scatter(x1,y1,'r','filled');hold on;
s2=scatter(x2,y2,'b','filled');hold on;
l=plot(xx,yfit,'Color','k','linewidth',2);hold on;
l=plot([0:0.1:1],polyval(p,[0:0.1:1]),'k--','linewidth',2);hold on;
legend([s1,s2,l],{'learner','non-learner',...
    ['linear fit: y=' num2str(p(1)) '*x+' num2str(p(2))]});
text(0.75,0.8,['R^2=' num2str(rsq)]);
axis equal;
xlim([0 1]);ylim([0 1]);
xlabel('Hab. shaking performance');
ylabel('Tst. learning performance');
set(gca,'fontsize',12);