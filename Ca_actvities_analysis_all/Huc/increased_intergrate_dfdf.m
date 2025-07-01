function [p,q,ind_increased]=increased_intergrate_dfdf(activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,env,outputpath)
global fs
global frame
global frameb
global trial
global re_startpoint
global startpoint

figure,
x = [1 2];
plot(x,area.CS_hab_tst,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
set(gca,'FontSize',20,'xtick',x,'xticklabel',{'Hab';'Test'},'TickLength',[0.01 0.01]);box off;
y_min = min(area.CS_hab_tst(:))-0.2*abs(min(area.CS_hab_tst(:)));
y_max = max(area.CS_hab_tst(:))+0.2*abs(min(area.CS_hab_tst(:)));
xlim([0.6 2.4]);ylim([y_min y_max]);
title('Area of CS','FontSize',20);
ylabel('Integrated △F/F0');

figure,
x = (1:1:size(area.CS_acq_block,1));
plot(x,area.CS_acq_block,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
xlim([0 size(area.CS_acq_block,1)]);
set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
title('Area of CS-US in each block','FontSize',20);
ylabel('Integrated △F/F0');xlabel('Blocks');

%page_test and 二项式拟合
permut=factorial(trial.acq_block_num);
p_page_acq=zeros(size(activities_preCS_dfdf_aftcorrect,1),1);
P_page_acq=zeros(size(activities_preCS_dfdf_aftcorrect,1),permut);
p_page_hab_tst=zeros(size(activities_preCS_dfdf_aftcorrect,1),1);
P_page_hab_tst=zeros(size(activities_preCS_dfdf_aftcorrect,1),permut);
for ii=1:size(activities_preCS_dfdf_aftcorrect,2)
    p.polyfit_acq(:,ii)=polyfit(1:trial.acq_block_num,area.CS_acq_block(:,ii)',1);
    p.polyfit_hab_tst(:,ii)=polyfit(1:2,area.CS_hab_tst(:,ii)',1);
    [p_page_acq(ii,1),P_page_acq(ii,:)] = mcpage(area.CS_acq_block(:,ii)',permut);
    p.pagetest_acq(:,ii)=mean(P_page_acq(ii,:)>=p_page_acq(ii));
    [p_page_hab_tst(ii),P_page_hab_tst(ii,:)] = mcpage(area.CS_hab_tst(:,ii)',permut);
    p.pagetest_hab_tst(:,ii)=mean(P_page_hab_tst(ii,:)>=p_page_hab_tst(ii));
end
% figure,histogram(P_page_acq(ind(1105),:),100);p.pagetest_acq(:,ind(1105))
% p_page_acq(ind(1105),1)
%adjust p cuase multi-test
q=p;
CL=0.05;
length(find(p.pagetest_acq<=CL))
sum(find(p.pagetest_hab_tst<=CL))

[q.pagetest_acq,fdr]=P_adj(p.pagetest_acq,'fdr');%'Storey','BHFDR';fdr矫正;Bonferroni校正
[q.pagetest_hab_tst,fdr]=P_adj(p.pagetest_hab_tst,'fdr');

length(find(q.pagetest_acq<=CL))
sum(find(q.pagetest_hab_tst<=CL))

%画图
ind_increased.pagetst_acq_nonadj=find(p.pagetest_acq<=CL);
plot_ind_mapback(p.polyfit_acq,ind_increased.pagetst_acq_nonadj,env,outputpath,['Acq_pagetest_p小于' num2str(CL)]) 

ind_increased.pagetst_acq=find(q.pagetest_acq<=CL);
plot_ind_mapback(p.polyfit_acq,ind_increased.pagetst_acq,env,outputpath,['Acq_pagetest_q小于' num2str(CL)]) 

ind_increased.pagetst_hab_tst=find(q.pagetest_hab_tst<=CL);
plot_ind_mapback(p.polyfit_acq,ind_increased.pagetst_hab_tst,env,outputpath,['Hab_tst_pagetest_q小于' num2str(CL)]) 


[~,thr,ind_increased.polyfit_acq]=plot_distribution(p.polyfit_acq,outputpath,['Acq_polyfit_k大于等于_CL' num2str(CL)],env,CL);
[~,thr,ind_increased.polyfit_hab_tst]=plot_distribution(p.polyfit_hab_tst,outputpath,['Hab_tst_k大于等于_CL' num2str(CL)],env,CL);


inter=intersect(ind_increased.pagetst_acq,ind_increased.polyfit_acq);
name=inter(zz(find(zz(:,1)~=0),1));%
plot_ind_mapback(p.polyfit_acq,name,env,outputpath,['Acq_intersect_pagetst_polyfit' num2str(CL)]); 
a=mean(area.CS_acq_block(:,name),2);
b=std(area.CS_acq_block(:,name),[],2);
a1=mean(area.CS_acq_block(:,setdiff(setdiff(inter,inter(zz(find(zz(:,2)~=0),2))),name)),2);%1:size(activities_preCS_dfdf_aftcorrect,2)
b1=std(area.CS_acq_block(:,setdiff(setdiff(inter,inter(zz(find(zz(:,2)~=0),2))),name)),[],2);

ind_increased_label=fieldnames(ind_increased);%{'ind_pagetst_acq','ind_pagetst_hab_tst','ind_polyfit_acq','ind_polyfit_hab_tst'};
%ind={ind_pagetst_acq,ind_pagetst_hab_tst,ind_polyfit_acq,ind_polyfit_hab_tst};
for kk=1:max(size(fieldnames(ind_increased)))
    path=checkpath([outputpath '\' ind_increased_label{kk}]);
    num=getfield(ind_increased,ind_increased_label{kk});
    if isempty(num)
        continue;warning([sheet{kk} 'is empty']);
    end
    checkpath([path '\heatmap']);
    checkpath([path '\raw trace']);
    acti=activities_preCS_dfdf_aftcorrect(:,num);
    a_event={};
    for ii=1:size(num,2)
        a_event{ii,1}=activities_event_preCS{num(ii)};
    end
%     for ii=1:1:size(num,2)
%         event=activities_event_preCS{num(ii),1};onset=[];
%         h1=figure;
%         plot_heatmap(acti(:,ii),event,fs,frame,trial,frameb,re_startpoint);
%         title([num2str(env.supervoxel(num(ii),3)) 'layer-' num2str(num(ii))]);xlabel('time(s)');ylabel('trial num.');colorbar;
%         saveas(h1,[path '\heatmap\' num2str(env.supervoxel(num(ii),3)) 'layer-' num2str(num(ii)) '.tif']);
%         close(h1);
%         
%         %hab;acq;test的raw trace
%         h2=figure;title([num2str(env.supervoxel(num(ii),3)) '-layer' num2str(num(ii)) ' Hab']);
%         plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,1);
%         saveas(h2,[path '\raw trace\' num2str(env.supervoxel(num(ii),3)) 'layer-' num2str(num(ii)) '.tif']);
%         close(h2);
%         
%         %spon bef and aft
%         if trial.spon_bef(1)~=0 && trial.spon_aft(1)~=0
%             checkpath([outpath '\raw trace-sp\']);
%             h3=figure;title([num2str(env.supervoxel(num(ii),3)) '-layer' num2str(num(ii))]);
%             plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,2);
%             saveas(h3,[path '\raw trace-sp\' num2str(env.supervoxel(num(ii),3)) 'layer-' num2str(num(ii)) '.tif']);
%             close(h3);
%         end
%     end
    y=mean(activities_preCS_dfdf_aftcorrect(:,num),2);event=getevent(y,trial,frame);
    %均值
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    plot_rawtrace(y,[1 frame.per_cycle*trial.total],[],fs,frame,trial,startpoint);legend('off');
    saveas(h,[path '\' 'mean' '.tif']);
    savefig(h,[path '\' 'mean']);
    close(h);
    %mean-heatmap
    h1=figure;
    plot_heatmap(y,event,fs,frame,trial,frameb,re_startpoint);
    title([ind_increased_label{kk} '-mean']);xlabel('time(s)');ylabel('trial num.');colorbar;
    saveas(h1,[path '\' '-mean-heatmap' '.tif']);
    savefig(h1,[path '\' '-mean-heatmap']);
    close(h1);
    %hab;acq;test的raw trace
    h2=figure;title([ind_increased_label{kk} '-mean' ' Hab']);
    plot_rawtrace_trials(y,event,fs,frame,trial,startpoint,1);
    saveas(h2,[path '\' 'mean-raw trace' '.tif']);
    savefig(h2,[path '\' 'mean-raw trace']);
    close(h2);
    %onset随时间变化
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    [summ1,summ2]=onset_errorbar(a_event,fs,trial,frame,re_startpoint);
    activities_event_acq_errorbar{kk,1}=[summ1.mx;summ1.m;summ1.sem]';
    activities_event_acq_errorbar{kk,2}=ind_increased_label{kk};
    weight=[[trial.acq(1):-1:1]/trial.acq(1)];
    activities_event_acq_errorbar{kk,3}=[summ2.ind;summ2.m;summ2.ind.*weight]';
    saveas(h,[path '\' 'errorbar onset' '.tif']);
    savefig(h,[path '\' 'errorbar onset']);
    close(h);
    %mean_acq
    hh1=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.4,0.7]) ;
    plot_meantrace_acq(y,fs,trial,frame);
    saveas(hh1,[path '\'  'meantrace_acq' '.tif']);
    savefig(hh1,[path '\'  'meantrace_acq']);
    close(hh1);
    %mean_hab_test
    hh2=figure;
    plot_mean_hab_test(y,fs,trial,frame,frameb,startpoint);
    saveas(hh2,[ path '\'  'trace_hab_test' '.tif']);
    savefig(hh2,[path '\'  'trace_hab_test']);
    close(hh2);
    
    %all traces sepplot
    ind=1:max(size(num));
    if ~isempty(ind)
        str=strcat({char(string(num(ind)'))});%,{char(char('-')*ones(size(num(ind)')))},{char(string(num(ind)'))}
        ind_frames=1:trial.total*frame.per_cycle;
        h=figure;
        set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.88]) ;
        plot_sepplot(ind_frames,acti(:,ind),str,frame,trial);hold on;
        saveas(h,[ path '\'  'all traces sepplot' '.tif']);
        savefig(h,[ path '\'  'all traces sepplot']);close(h);
        %figure,Correlational_analysis(acti(:,ind));
        h=figure;subplot(1,2,1);
        Correlational_analysis(acti((trial.hab(2)-1)*frame.per_cycle+1:trial.hab(3)*frame.per_cycle,ind));
        subplot(1,2,2);
        Correlational_analysis(acti((trial.test(2)-1)*frame.per_cycle+1:trial.test(3)*frame.per_cycle,ind));
        saveas(h,[ path '\'  'all traces corr_hab_test' '.tif']);
        savefig(h,[ path '\'  'all traces corr_hab_test']);
        close(h);
    end
end
% h=figure;
% [m,thr]=plot_distribution(p.polyfit_hab_tst,[],'<');
% saveas(h,[outputpath '\'  'slope_distribution_hab_tst_less' '.tif']);
% savefig(h,[outputpath '\'  'slope_distribution_hab_tst_less']);
% vmap = mapback(abs(p.hab_tst(1,ind)), env.supervoxel,[env.height env.width env.depth], ind);
% seqwrite(vmap, checkpath([outputpath '/Acq_k小于等于' num2str(thr)]));