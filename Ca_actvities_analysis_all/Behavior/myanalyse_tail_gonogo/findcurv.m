function [curv,curv1]=findcurv(rawdata)
%≤Â÷µor≤ª≤Â÷µ
% n=200;
% for i=1:size(rawdata,3)
% rawdata(1:n,1,i)=interp1(rawdata(:,2,i),rawdata(:,1,i),linspace(min(rawdata(:,1,i)),max(rawdata(:,1,i)),n),'pchip');
% rawdata(1:n,2,i)=interp1(rawdata(:,1,i),rawdata(:,2,i),linspace(min(rawdata(:,2,i)),max(rawdata(:,2,i)),n),'pchip');
% end
%pointi-pionti+1  i+1curv
vec=rawdata(1:end-1,:,:)-rawdata(2:end,:,:);
curv(:,:)=atan2d(vec(:,2,:),vec(:,1,:));
curv1(:,:)=atand(vec(:,2,:)./vec(:,1,:));
end