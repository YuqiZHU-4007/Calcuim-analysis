 %%%%%%%%%%%%% pause��ťcallback���� %%%%%%%%%%%%%%%%%
    function callback_pause(hobj,dataevent,stopexampling)
         hfig=findobj('tag','hmainfigure');
         hrunpush=findobj('tag','run_push');
         currenttest=get(findobj(hfig,'tag','curitem_popup'),'value');
         set(hrunpush,'string','Run');
          if ~isempty( timerfind('tag','testtimer'))%�ر�timer
              stop(timerfind('tag','testtimer'));
              delete(timerfind('tag','testtimer'));
           end
         if stopexampling==1 & currenttest~=1 %search cellʱ ��pause��ť���ᷢ�жϲɼ����ź�
           pause(0.3);
           crsIOWriteDigitalOut(0,CRS.DIG2);%�ѵ�3��Ÿ�λ���͵�ƽ
         end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%