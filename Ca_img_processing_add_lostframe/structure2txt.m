%%此函数用于把参数结构体存储为txt格式%%%%%%%%%%%%%%%%%%%
function structure2txt(struc,strucname,txtname)
%%struc is the structure.(its form is structure)
%%strucname is the name of structure that to be saved.(its form is txt)
%%txtname is the file name as which the structure will be saved.(its form
%%is txt)
paracell=fieldnames(struc);
fid=fopen(txtname,'wt');
for i=1:length(paracell)
    currentvalue=getfield(struc,cell2mat(paracell(i)));
    if length(currentvalue)==1
        fprintf(fid,[strucname '.' cell2mat(paracell(i)) '=%s;\n'],num2str(currentvalue));
     else
        fprintf(fid,[strucname '.' cell2mat(paracell(i)) '=[%s];\n'],num2str(currentvalue,'%6.2f '));
    end
end