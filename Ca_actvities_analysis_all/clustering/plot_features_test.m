function plot_features_test(x,cIX,gIX,ref_win,win,fs)

for i = 1:length(unique(gIX))
    
    clust = find(gIX==i);
    xx=mean(x(:,cIX(clust)),2);
    para=getpara_event_20190521(xx,ref_win,win,fs);
    name=fieldnames(para);
    txt_dis=string;
    for ii=1:length(name)
        a = getfield(para,name{ii});
        txt_dis=[txt_dis;name{ii},': ',num2str(a)];
    end
    figure,
    txt = uicontrol('Style','text',...
        'Position',[20,1,200,400],...
        'fontsize',15,...
        'BackgroundColor','w',...
        'horizontalalignment','left',...
        'string',txt_dis);
end