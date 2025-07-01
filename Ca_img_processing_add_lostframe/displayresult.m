function displayresult(a)
str={};
hfc=findobj('tag','hmainfigure');
hdis=findobj(hfc,'tag','dispanel');
str=get(findobj(hdis,'tag','dis_tag'),'string');
if iscell(str)
%    disstr=[str;a];
   % set(findobj(hdis,'tag','dis_tag'),'string',[str;a]);
    switch class(a)
        case 'char'
            set(findobj(hdis,'tag','dis_tag'),'string',([str;a]));
            %set(findobj(hdis,'tag','dis_tag'),'string',strcat(str, {', '}, a));
        otherwise
            set(findobj(hdis,'tag','dis_tag'),'string',([str;struct2str(a)]));
    end
else
    switch class(a)
        case 'char'
            set(findobj(hdis,'tag','dis_tag'),'string',cellstr([str;a]));
            %set(findobj(hdis,'tag','dis_tag'),'string',strcat(str, {', '}, a));
        otherwise
            set(findobj(hdis,'tag','dis_tag'),'string',cellstr([str;struct2str(a)]));
    end
end
end