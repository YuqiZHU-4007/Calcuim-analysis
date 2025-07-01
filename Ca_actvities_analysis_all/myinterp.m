function yy=myinterp(x,y,q,method)
%20180926修改版,因为可能动的那点才是最高值，如果线性插值，会让最高点后移，导致判断不到反应

%x,y,q:1-d
yy=y;
for i=1:length(q)
    a=x(find(x<q(i)));a=max(a);
    b=x(find(x>q(i)));b=min(b);
    
    switch method
        case 'linear'%'linear'
            yy(q(i))=(q(i)-a)*(y(b)-y(a))/(b-a)+y(a) ;
        case 'max' %插成最大值
            yy(q(i))=max(y(b),y(a)) ;
    end
end