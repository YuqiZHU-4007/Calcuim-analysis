 
function [gIX,C,numK,clrmap]=find_best_k_in_range(M,range)
% gIX = getappdata(hfig,'gIX');
% M = getappdata(hfig,'M');

if ~isempty(range),
    if length(range)~=1
    for i = 1:length(range),
        disp(['k-means k=' num2str(range(i))]);
        [gIX,~,sum ]= kmeans(M,range(i),'distance','correlation','Replicates',100,'MaxIter',1000,'Display','final');
        silh = silhouette(M,gIX,'correlation');
        Silh_mean(i) = mean(silh);
        %sumD(i)=sum;
    end
    [~,ix] = max(Silh_mean);
    numK = range(ix);
    
    figure;hold on
    plot(range,Silh_mean,'o-','color',[0.5 0.5 0.5]);
    
    if length(range)~=1,
        plot([numK,numK],[min(Silh_mean),max(Silh_mean)],'r--');
        xlim([range(1),range(end)]);
        ylim([min(Silh_mean),max(Silh_mean)])
    end
    xlabel('k = # of clusters')
    ylabel('silhouette mean')
    box on;
    else
        numK=range;
    end
%     figure;
%     plot(range,sumD,'o-','color',[0.5 0.5 0.5]);
    %% perform regular k-means
    disp(['k-means k=' num2str(numK)]);

    %rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK < 10^7 && numK~=1,
        disp('Replicates = 5');
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        [gIX,C] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
    end
    if numK>1,
        disp('HierClus_Direct');
        gIX = HierClus_Direct(C,gIX);
    end
    if ~isempty(numK)
        clrmap = GetColormap('hsv_new',numK);
    else
        clrmap = [];
    end
end
end