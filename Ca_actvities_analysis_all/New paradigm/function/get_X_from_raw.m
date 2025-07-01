function [A,B]=get_X_from_raw(n,raw_data,rownum_eachfish)
%A:wirte to excel,graphpad
 A=[];B=[];
if length(n)==1
    for jj=1:ceil(size(raw_data,2)/rownum_eachfish)
        A=[A,raw_data(:,n+rownum_eachfish*(jj-1)),nan(size(raw_data(:,n+rownum_eachfish*(jj-1)))),nan(size(raw_data(:,n+rownum_eachfish*(jj-1))))];
        B=[B raw_data(:,n+rownum_eachfish*(jj-1))];
    end
elseif length(n)==2
    for jj=1:ceil(size(raw_data,2)/rownum_eachfish)
        A=[A,raw_data(:,n+rownum_eachfish*(jj-1)),nan(size(raw_data(:,n+rownum_eachfish*(jj-1)),1),1)];
        B=[B raw_data(:,n+rownum_eachfish*(jj-1))];
    end
end
if  sum(n==[8:9])==2 || sum(n==3)==1
    A(3:end,:)=[];
    B(3:end,:)=[];
end