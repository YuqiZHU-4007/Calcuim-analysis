function [In,On]=count_spatial_location_density(loc,all_loc, radius)
% %2d
% %radius£º°ë¾¶
% for ii=1:length(loc)
%     r = linspace(0,2*pi,100);
%     xv = radius*cos(r)'+loc(ii,1);
%     yv = radius*sin(r)'+loc(ii,2);
%     xq=all_loc(:,1);yq=all_loc(:,2);
%     [in,on] = inpolygon(xq,yq,xv,yv);
%         numel(xq(in))
%         numel(xq(on))
%         figure,plot(xv,yv) % polygon
%         axis equal;hold on
%         plot(xq(in),yq(in),'r+') % points inside
%         xlim([0 2048]);ylim([0 2048]);
%         plot(xq(~in),yq(~in),'k.') % points outside
%         hold off
% end

%3d
%radius£º°ë¾¶
In=zeros(size(loc,1),1);
On=zeros(size(loc,1),1);
for ii=1:size(loc,1)
    r = linspace(0,2*pi,100);
    xv = radius*cos(r)'+loc(ii,1);
    yv = radius*sin(r)'+loc(ii,2);
    ind=find(all_loc(:,3)==loc(ii,3));
    xq=all_loc(ind,1);yq=all_loc(ind,2);
    [in,on] = inpolygon(xq,yq,xv,yv);
    In(ii)=    numel(xq(in));
    On(ii)=    numel(xq(on));
    %         figure,plot(xv,yv) % polygon
    %         axis equal;hold on
    %         plot(xq(in),yq(in),'r+') % points inside
    % %         plot(xq(~in),yq(~in),'bo') % points outside
    %         hold off
end

map = colormap(parula(max(In)));
ind = In;fix((In-min(In))/(max(In)-min(In))*(size(map,1)-1))+1;
%%3d
% kk=1;
% figure,
% for ii=unique(loc(:,3))'
%     subplot(4,6,kk)
%     if any(loc(:,3)==ii)
%         ind_z=find(loc(:,3)==ii);
%         gsp(loc(ind_z,:),ind(ind_z),map);
%     end
%     kk=kk+1;
% end
%%2d
%figure,
gsp(loc,ind,map);

    function h=gsp(loc,ind,map)
        h = [];
        for k=1:size(map,1)
            if any(ind==k)
                h(end+1) = line('Xdata',loc(ind==k,1),'Ydata',loc(ind==k,2), 'Zdata',loc(ind==k,3),...
                    'LineStyle','none','Color',map(k,:), ...
                    'Marker','.','MarkerSize',3);
            end
        end
        if nargout==1
            varargout{1} = h;
        end
        %colorbar('Ticks',[0,0.5,1],'TickLabels',[0 max(In)/2 max(In)]);
        axis equal;xlim([0 2048*0.66]);ylim([0 2048*0.66]);%set(gca,'visible','off')
    end
end
