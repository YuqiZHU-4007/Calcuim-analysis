function yy=myinterp(x,y,q,method)
%20180926�޸İ�,��Ϊ���ܶ����ǵ�������ֵ��������Բ�ֵ��������ߵ���ƣ������жϲ�����Ӧ

%x,y,q:1-d
yy=y;
for i=1:length(q)
    a=x(find(x<q(i)));a=max(a);
    b=x(find(x>q(i)));b=min(b);
    
    switch method
        case 'linear'%'linear'
            yy(q(i))=(q(i)-a)*(y(b)-y(a))/(b-a)+y(a) ;
        case 'max' %������ֵ
            yy(q(i))=max(y(b),y(a)) ;
    end
end