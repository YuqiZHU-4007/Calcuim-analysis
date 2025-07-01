function callback_stop
global key
%%stop 按钮callback程序
pause(5)
hfig=findobj('tag','hmainfigure');
%%关闭呈现刺激的timer
        if ~isempty( timerfind('tag','testtimer'))%关闭timer
           stop(timerfind('tag','testtimer'));
           delete(timerfind('tag','testtimer'))
        end
%        set(findobj(hfig,'string','d'),'enable','on');
key=0;
end