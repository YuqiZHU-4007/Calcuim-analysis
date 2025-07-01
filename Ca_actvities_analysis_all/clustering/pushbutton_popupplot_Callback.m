function [h]=pushbutton_popupplot_Callback(M,cIX,gIX,clrmap,env,fs,stimCS,stimUS,isPlotAnatomyOnly,i_fish,isPlotLines,isPlotBehavior,isCentroid,isPlotRegWithTS)
% very similar as function RefreshFigure(hfig)
isPopout = 1; % no down-sampling in plots
isLeftPlotOnly = 1; % manual...

if ~isPlotAnatomyOnly
    h=figure('Position',[50,50,1200,600],...  % [50,50,1600,800],
        'color',[1 1 1],...
        'Name',['Fish#' num2str(i_fish)]);
    h1 = axes('Position',[0.05, 0.03, 0.6, 0.94]); % left ~subplot
    h2 = axes('Position',[0.66, 0.03, 0.32, 0.94]); % right ~subplot
    
    % left subplot
    axes(h1);
    opts = [];
    opts.h_ax = h1;
    opts.i_fish=i_fish;
    opts.isPopout = isPopout;
    opts.isCentroid = isCentroid;
    opts.isPlotLines = isPlotLines;
    opts.isPlotBehavior = isPlotBehavior;
    opts.isPlotRegWithTS = isPlotRegWithTS;
    opts.numK=length(unique(gIX));
    opts.stimCS=stimCS;
    opts.stimUS=stimUS;
    DrawTimeSeries_zyq_20190530(M(cIX,:),cIX,gIX,clrmap,fs,opts)
    
    %     % right subplot
    axes(h2);
    value=repclrmap(clrmap,gIX);
    showspv= mapback(value, env.supervoxel,[env.height env.width 1],cIX);
    imshow(showspv)
%     circlemaskIX = MakeCircularMask(5,[env.width,env.height,3]); 
%     % get colormap
%     numK = opts.numK;
%     if isempty(numK),
%         numK = max(gIX);
%     else
%         numK = double(max(numK,max(gIX)));% sadly this is not always true
%     end
%    
%     % set transparancy
%     alpha_max = 0.4;%0.3-min(length(cIX)/1000/100,0.1);
%     clr_alpha = ones(size(cIX))*alpha_max;
%     vol(:,:,1)=env.vol(:,:,1);vol(:,:,2)=env.vol(:,:,1);vol(:,:,3)=env.vol(:,:,1);
%     anat_YX=DrawMasksInRGB(vol,env.supervoxel(cIX,2:1),circlemaskIX,clrmap,gIX,clr_alpha);
%     imshow(anat_YX,[min(anat_YX(:)) max(anat_YX(:))]); % image(tot_image);
%     axis image;axis off
%     set(gcf,'Color','w');
    %     DrawCellsOnAnatProj(hfig,isRefAnat,isPopout);
elseif isLeftPlotOnly
    h=figure('Position',[412,225,414,175],... % [50,50,1600,800],[50,50,600,950],
        'color',[1 1 1],...
        'Name',['Fish#' num2str(i_fish)]);
    h1 = axes('Position',[0.1, 0.03, 0.85, 0.94]); % left ~subplot
    %     h2 = axes('Position',[0.46, 0.03, 0.5, 0.94]); % right ~subplot
    
    % left subplot
    axes(h1);
    opts = [];
    opts.h_ax = h1;
    opts.i_fish=i_fish;
    opts.isPopout = isPopout;
    opts.isCentroid = isCentroid;
    opts.isPlotLines = isPlotLines;
    opts.isPlotBehavior = isPlotBehavior;
    opts.isPlotRegWithTS = isPlotRegWithTS;
    opts.numK=length(unique(gIX));
    opts.stimCS=stimCS;
    opts.stimUS=stimUS;
    DrawTimeSeries_zyq_20190530(M(cIX,:),cIX,gIX,clrmap,fs,opts)
else
    %     % right subplot
    %     I = LoadCurrentFishForAnatPlot(hfig);
    %     DrawCellsOnAnat(I);
    % %     AddColorbarToAnat(I.clrmap);%,cmin,cmax)
end
end
