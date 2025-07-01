function name_fishnumber=get_fishnumber(raw_data,rownum_eachfish)
kk=1;
for jj=1:ceil(size(raw_data,2)/rownum_eachfish)
    %             name{1,1+3*(jj-1)}={raw_data{1,1+rownum_eachfish*(jj-1)}};
    %             name{1,1+3*(jj-1)+1}=nan;
    %             name{1,1+3*(jj-1)+2}=nan;
    if isnan(raw_data{1,1+rownum_eachfish*(jj-1)})
        continue;
    else
        name_fishnumber{1,kk}=raw_data{1,1+rownum_eachfish*(jj-1)};
        kk=kk+1;
    end
end
