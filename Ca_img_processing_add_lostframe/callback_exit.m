function callback_exit
%%�˳���������
global key
key=0;
close(gcf);
if ~isempty( timerfind('tag','testtimer'))%�ر�timer
    stop(timerfind('tag','testtimer'));
    delete(timerfind('tag','testtimer'));
end
button = questdlg('Ready to quit?','Exit Dialog','Yes','No','No');
switch button
    case 'Yes',
        disp('Exiting MATLAB')
        % Save variables to matlab.mat
        %save(savepath)
        exit 
    case 'No',
        quit cancel;
end
end