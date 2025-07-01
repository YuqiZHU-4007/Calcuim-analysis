a=3;
index=ind_increased.intersect_acq_pagetst_polyfit_increased(a);
figure,
imshow(env.vol(:,:,env.supervoxel(index,3)),[0 500]);
hold on;scatter(env.supervoxel(index,1),env.supervoxel(index,2),env.supervoxel(index,5),'r','filled');
%scatter(2048,1,env.supervoxel(index,4),'r','filled');
figure,
subplot(2,1,1),plot([activities(index,1) activities(index,:)]);xlim([1 size(activities,2)+1]);ylabel('F');set(gca,'fontsize',20);
subplot(2,1,2),plot([activities_preCS_dfdf_aftcorrect(:,index)]);xlim([1 size(activities,2)+1]);ylabel('dfdf');set(gca,'fontsize',20);