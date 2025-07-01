function plot_ind_mapback(fit,ind,env,outputpath,tag) 
if ~isempty(ind)
    %figure,plot(area.CS_acq_block(:,ind));
    %fit=p.polyfit_acq;
    h=figure;
    scatter(fit(1,:),fit(2,:),[],'b');hold on;
    pp=[];
    pp(1,:)=fit(1,ind);pp(2,:)=fit(2,ind);
    scatter(pp(1,:),pp(2,:),[],'r');hold on;
    xlabel('Slope');ylabel('Intercept');
    set(gca,'fontsize',20);
    title(['slope_distribution_' tag],'Interpreter','none');
    saveas(h,[outputpath '\slope_distribution_' tag '.tif']);
    savefig(h,[outputpath '\slope_distribution_' tag '.fig']);
    vmap = mapback(fit(1,ind), env.supervoxel,[env.height env.width env.depth], ind);
    seqwrite(vmap, checkpath([outputpath '\' tag]));
else
    disp('can not find singnificant acq increased neuron');
end