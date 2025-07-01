function [h,gIX2, numU] = hierplot_zyq_20190530(cIX,gIX,M)

[h,gIX2, numU] = HierClus(M,gIX,'isplotfig');
if ~isequal(gIX,gIX2),
    disp('gIX2 is not equal to gIX');
%     UpdateIndices(hfig,cIX,gIX2,numU);
%     RefreshFigure(hfig);
end
end

function [h,gIX, numU] = HierClus(M,gIX,isplotfig) %#ok<INUSD>
[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid_Direct(gIX,M);
D = pdist(C,'correlation');
if size(C,1)>1,
    tree = linkage(C,'average','correlation');
    leafOrder = optimalleaforder(tree,D);
    
    if numU>1,
        if exist('isplotfig','var'),
            h=figure('Position',[100 100 600 600]);
            %             subplot(1,3,1);
            %             CORR = corr(C');
            %             CorrPlot(CORR);
            %
            %             subplot(1,3,2);
            %dendrogram(tree,0,'orientation','right');
            dendrogram(tree,numU,'orientation','right','reorder',leafOrder);
            set(gca,'YDir','reverse');
%             set(gca,'XTick',[]);
            
            %             subplot(1,3,3);
            %             C2 = C(leafOrder,:);
            %             CORR2 = corr(C2');
            %             CorrPlot(CORR2);
        end
        % sort for uniform colorscale
        temp = zeros(size(gIX));
        for i = 1:numU,
            temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
        end
        gIX = temp;
    end
end
end
