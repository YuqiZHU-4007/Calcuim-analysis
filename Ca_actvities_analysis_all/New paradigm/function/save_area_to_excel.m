function  save_area_to_excel(S,filepath,batchi,fishname)

fields = fieldnames(S) ;

for ii=1:length(fields)
    sheetname=getfield(S,'name');
    if strcmpi(fields{ii},'name')
        continue;
    else
        path_this_field=checkpath([filepath '\' fields{ii}]);
        A=getfield(S,fields{ii});
        XLrange=char(65+(batchi-1)*(size(A,2))/length(sheetname))
        write_A_to_xlsfile(A,sheetname,path_this_field,batchi,fishname)
    end
end