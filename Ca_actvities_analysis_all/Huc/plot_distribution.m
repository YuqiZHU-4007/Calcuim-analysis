function [m,thr,ind]=plot_distribution(a,outputpath,tag,env,thr,thr_type)

m=mean(a(1,:));sd=std(a(1,:));
if nargin<5 || isempty(thr)
    thr=m+2*sd;
    thr_type='>';
elseif nargin<6 || isempty(thr)
    thr=m+2*sd;
    thr_type='>';
end
switch thr
    case 0.01
        thr=m+3*sd;
    case 0.05
        thr=m+2*sd;
end
h1=figure;
subplot(4,2,1);
scatter(a(1,:),ones(size(a(1,:))),[],[0.8,0.08,0.18]);hold on;
set(gca,'fontsize',20,'xcolor','w','ycolor','w');
subplot(4,2,[3 5 7]);
h=histogram(a(1,:),160,'FaceColor',[0.8,0.08,0.18]);hold on;
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,h.Values,'r','linewidth',2);hold on;
%h2=histogram([p.acq(1,find(p.acq(1,:)<=0)) abs(p.acq(1,find(p.acq(1,:)<=0)))],h.BinEdges,'FaceColor',[0.0,0.45,0.74]);
%h=histogram(p.acq(1,:),160,'FaceColor',[0.8,0.08,0.18]);hold on;
line([m m],[0 max(h.Values)],'color','k','linewidth',1.5,'linestyle','--');hold on;
line([thr thr],[0 max(h.Values)],'color','b','linewidth',1.5,'linestyle','-');hold on;
%subplot(5,1,5);bar((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,h.Values./h2.Values);
set(gca,'fontsize',20);box off

switch thr_type
    case '>'
        ind=find(a(1,:)>=thr);
    case '<'
        ind=find(a(1,:)<thr);
end
subplot(4,2,[2 4 6 8]);
scatter(a(1,:),a(2,:),[],'b');hold on;
pp=[];
pp(1,:)=a(1,ind);pp(2,:)=a(2,ind);
scatter(pp(1,:),pp(2,:),[],'r');hold on;
xlabel('Slope');ylabel('Intercept');
set(gca,'fontsize',20);
title(['slope_distribution_' tag],'Interpreter','none');

saveas(h1,[outputpath '\slope_distribution_' tag '.tif']);%'\slope_distribution_acq_ind_polyfit_acq' '.tif'
savefig(h1,[outputpath '\slope_distribution_' tag '.fig']);
vmap = mapback(a(1,ind), env.supervoxel,[env.height env.width env.depth], ind);
seqwrite(vmap, checkpath([outputpath '\' tag]));%'/Acq_polyfit_k大于等于2sd-' num2str(thr) '_CL' num2str(CL)