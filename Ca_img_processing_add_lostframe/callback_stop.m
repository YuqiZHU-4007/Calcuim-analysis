function callback_stop
global key
%%stop ��ťcallback����
pause(5)
hfig=findobj('tag','hmainfigure');
%%�رճ��ִ̼���timer
        if ~isempty( timerfind('tag','testtimer'))%�ر�timer
           stop(timerfind('tag','testtimer'));
           delete(timerfind('tag','testtimer'))
        end
%        set(findobj(hfig,'string','d'),'enable','on');
key=0;
end