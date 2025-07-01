function fig=show_spv_GUI(image)

K=size(image);K=K(end);
fig = figure('Visible','off');
set(gcf,'Position',2*[100,30,420,450]);
set(gcf,'PaperPosition',2*[300,300,960,480]);
% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',K,'Value',1,'SliderStep',[1/(K-1) 1],...
    'Position', [150 20 420 20],...
    'Callback', @surfzlim);
% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 150 20],...
    'String','Slice');
% Make figure visble after adding all components
fig.Visible = 'on';

dim=length(size(image));
if dim==3
    a=image(:,:,1);
    imshow(image(:,:,1),[min(a(:)),max(a(:))]);hold on;
elseif dim==4
    a=image(:,:,:,1);
    if min(a(:))>= max(a(:))
        imshow(image(:,:,:,1),[0,0.5]);hold on;
    else
        imshow(image(:,:,:,1),[min(a(:)),max(a(:))]);hold on;
    end
end

    function surfzlim(source,callbackdata)
        i = source.Value;
        set(txt,'String',['Slice ' num2str(i)]);
        cla;
        if dim==3
            a=image(:,:,round(i));
            imshow(image(:,:,round(i)),floor([min(a(:)),max(a(:))]));hold on;
        elseif dim==4
            a=image(:,:,:,round(i));
            imshow(image(:,:,:,round(i)),([min(a(:)),max(a(:))]));hold on;
        end
    end
hold off
end
