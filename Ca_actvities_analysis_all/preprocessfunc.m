function x=preprocessfunc(x,type)
s=size(x);
x=reshape(x,1,[]);
switch type
    case 1
        x = scaleMatrixToRange(x, -1, 1);     
    case 2
        x=normalize(x,'zscore');
    case 3
        x=normalize(x,'range');       
end
%x=(x-min(x))./(max(x)-min(x));
x=reshape(x,s(1),s(2),[]);
x = smoothdata(x,2);
end
