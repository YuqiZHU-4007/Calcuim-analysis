function mapback_colorcode_responses(loc,ind)
% %2d
if ~isempty(loc)
    %3d
    map = colormap(hot(length(unique(ind))));
    gsp(loc,ind,map);
end
    function h=gsp(loc,ind,map)
        h = [];a=unique(ind);
        for k=1:length(a)
            if any(ind==a(k))
                h(end+1) = line('Xdata',loc(ind==a(k),1),'Ydata',loc(ind==a(k),2), 'Zdata',loc(ind==a(k),3),...
                    'LineStyle','none','Color',map(k,:), ...
                    'Marker','.','MarkerSize',3);
            end
        end
        if nargout==1
            varargout{1} = h;
        end
        colorbar('Ticks',[0,0.5,1],'TickLabels',[0 max(ind)/2 max(ind)]);
        axis equal;xlim([0 2048]);ylim([0 2048]);set(gca,'visible','off')
    end
end


