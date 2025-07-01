 %%%%%%%%%%%%% pause按钮callback函数 %%%%%%%%%%%%%%%%%
    function callback_pause(hobj,dataevent,stopexampling)
         hfig=findobj('tag','hmainfigure');
         hrunpush=findobj('tag','run_push');
         currenttest=get(findobj(hfig,'tag','curitem_popup'),'value');
         set(hrunpush,'string','Run');
          if ~isempty( timerfind('tag','testtimer'))%关闭timer
              stop(timerfind('tag','testtimer'));
              delete(timerfind('tag','testtimer'));
           end
         if stopexampling==1 & currenttest~=1 %search cell时 按pause按钮不会发中断采集的信号
           pause(0.3);
           crsIOWriteDigitalOut(0,CRS.DIG2);%把第3针脚复位到低电平
         end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%