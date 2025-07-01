load('W:\RNN\Post_20210709fish2.mat');
flow=nan(size(J,1),size(J,2),size(trace_data,2));
for ii=1:size(J,1)
    for jj=1:size(J,2)
        flow(ii,jj,:)=J(ii,jj)*trace_data(ii,:);
    end
end

regioni=3;
source= double(regions{regioni,2})+1;
weight=J0_r(1,source);
 flow_region=nan(length(regions),size(trace_data,2));
for regionj=1:length(regions)
    target=double(regions{regionj,2})+1;
    flow_region_i= squeeze(nanmean(flow(target,source,:),1));
    flow_region_j=[];
    for jj=1:length(source)
        flow_region_j(jj,:)=  flow_region_i(jj,:).*weight(jj);
    end
    flow_region(regionj,:)=nanmean(flow_region_j,1);
end

figure,
for regionj=1:length(regions)
subplot(5,5,regionj),plot(flow_region(regionj,:)*100000);title(regions{regionj,1});
ylim([-5 10])
end