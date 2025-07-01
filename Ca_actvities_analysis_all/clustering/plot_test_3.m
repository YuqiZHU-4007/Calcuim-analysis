function [h1,h2,ratio_fish]=plot_test_3(A,cIX,gIX,clrmap,index_all,envbatch,env_all,outputpath,fs,stimCS,stimUS)
ratio_fish=[];h1={};
for kk=unique(index_all)'
    [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX);
    ind_cls_in_this_fish=i_cIX;
%     %ind_cls_in_this_fish=intersect(CR_loc(cIX2),ind_event_in_this_fish);
%     ind_cls_in_this_fish=[];zz=1;
%     for ii=1:length(ind_event_in_this_fish)
%         iinddd=find(cIX==ind_event_in_this_fish(ii));
%         if ~isempty(iinddd)
%             ind_cls_in_this_fish(zz)=iinddd;
%             zz=zz+1;
%         end
%     end
    load(envbatch{kk});
    act_trace=A';
    for ii=unique(gIX(ind_cls_in_this_fish))'
        ratio_fish(kk,ii)=length(find(gIX(ind_cls_in_this_fish)==ii));
    end
    colorind=unique(gIX(ind_cls_in_this_fish));
%     [h1]=pushbutton_popupplot_Callback(act_trace,cIX(ind_cls_in_this_fish),gIX(ind_cls_in_this_fish),clrmap(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
%     saveas(h1,[outputpath '\fish' num2str(kk,'%02d') '_1'],'fig');
    h2=DrawTiledPics_zyq_20190530(cIX(ind_cls_in_this_fish),gIX(ind_cls_in_this_fish),[1:size(act_trace,1)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap);
    %saveas(h2,[outputpath '\fish' num2str(kk,'%02d') '_2'],'fig');
    %figure,scatter(env_all.supervoxel(cIX(ind_cls_in_this_fish),1),env_all.supervoxel(cIX(ind_cls_in_this_fish),2));
end

% for kk=unique(index_all)
%     [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX);
%     ind_cls_in_this_fish=i_cIX;
% %     %ind_cls_in_this_fish=intersect(CR_loc(cIX2),ind_event_in_this_fish);
% %     ind_cls_in_this_fish=[];zz=1;
% %     for ii=1:length(ind_event_in_this_fish)
% %         iinddd=find(cIX==ind_event_in_this_fish(ii));
% %         if ~isempty(iinddd)
% %             ind_cls_in_this_fish(zz)=iinddd;
% %             zz=zz+1;
% %         end
% %     end
%     load(envbatch{kk});
%     act_trace=act_all_acq';
%     for ii=unique(gIX(ind_cls_in_this_fish))'
%         ratio_fish(kk,ii)=length(find(gIX(ind_cls_in_this_fish)==ii));
%     end
%     colorind=unique(gIX(ind_cls_in_this_fish));
%     [h]=pushbutton_popupplot_Callback(act_trace,cIX(ind_cls_in_this_fish),gIX(ind_cls_in_this_fish),clrmap(colorind,:),env,fs.ca,stimCS,stimUS,1,kk,0,0,0,0);
%     saveas(h,[outputpath '\fish' num2str(kk,'%02d') '_1'],'fig');
%     h=DrawTiledPics_zyq_20190530(cIX(ind_cls_in_this_fish),gIX(ind_cls_in_this_fish),[1:size(act_trace,1)],[env_all.supervoxel(:,2) env_all.supervoxel(:,1) env_all.supervoxel(:,3)],env.vol,clrmap);
%     saveas(h,[outputpath '\fish' num2str(kk,'%02d') '_2'],'fig');
% end