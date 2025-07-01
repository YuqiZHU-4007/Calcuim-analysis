function [interA,interB,inter]=findintersect(A,B,bin)

if nargin<3
    bin=32;
end
win=ceil(bin/2);kk=1;
for ii=1:size(A,1)
    %b=B(find(B(:,1)==A(ii,1)),2);
    ind=find(B(:,2)<A(ii,2)+win & B(:,2)>=A(ii,2)-win & B(:,1)==A(ii,1));
    if ~isempty(ind) && length(ind)==1
        interA(kk,:)=A(ii,:);
        interB(kk,:)=B(ind,:);
        inter{kk,1}=A(ii,:);
        inter{kk,2}=B(ind,:);
        kk=kk+1;
    end
end




% figure,
% x1 = 0; y1 = 0; r1 = length(A)/(length(A)+length(B));
% x2 = 1; y2 = 0; r2 = length(B)/(length(A)+length(B));
% [xout,yout] = circcirc(x1,y1,r1,x2,y2,r2);
% f1 =@(x) sqrt( -(x-x1).^2 + (r1)^2 ) + y1;
% f2 =@(x) -sqrt( -(x-x1).^2 + (r1)^2 ) + y1;
% f3 =@(x) sqrt( -(x-x2).^2 + (r2)^2 ) + y2;
% f4 =@(x) -sqrt( -(x-x2).^2 + (r2)^2 ) + y2;
% x = xout(1):(xout(end)-xout(1))/50:xout(end);
% X1 = (x1-r1):(r1*2)/50:(x1+r1);
% X2 = (x2-r2):(r2*2)/50:(x2+r2);
% plot(X1,f1(X1),'r-',X1,f2(X1),'r-',X2,f3(X2),'b-',X2,f4(X2),'b-');hold on;
% fill([x x],[f1(x) f4(x)],'k','LineStyle','none')