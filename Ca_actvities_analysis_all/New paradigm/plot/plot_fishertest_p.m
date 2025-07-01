function [h,p,stats]=plot_fishertest_p(f,xx,yy,catege,higher_than_edge)
[h,p,stats] = fishertest(f);
if h==1
    %¾Ü¾øÁã¼ÙÉè£¬ÓĞÏÔÖøĞÔ²îÒì
    line(catege,[max([xx,yy]) max([xx,yy])]+higher_than_edge,'Color', 'k');hold on;
    scatter(catege(1),max([xx,yy])+higher_than_edge*2,15,'k','*');hold on;
else
    line(catege,[max([xx,yy]) max([xx,yy])]+higher_than_edge,'Color', 'k');hold on;
    text(catege(1),max([xx,yy])+higher_than_edge*2,15,'ns');hold on;
end