function [taildeg_everytrial,taildeg_everytrial1]=mov_detection_1(alltheta,bin1,bin2,trailnum)
alltheta(alltheta > 180) = alltheta(alltheta > 180) - 360;
alltheta(alltheta < -180) = alltheta(alltheta < -180) + 360;
for i=1:trailnum
    taildeg_everytrial(i,:)=alltheta(1+12000*(i-1):12000+12000*(i-1),:);
%     taildeg_sti(i,:)=alltheta(7201+12000*(i-1):7440+12000*(i-1),:);
%     taildeg_bese(i,:)=alltheta(1+12000*(i-1):sti_start+12000*(i-1),:);
    % a=find((taildeg_sti(i,:)<5 & 0<taildeg_sti(i,:))|(taildeg_sti(i,:)<0 & -5<taildeg_sti(i,:)));
    % if a~=0
    % a1=[a(2:end) a(end-1)];aa=a1-a;%%说明第二个时刻也在阈值内
    % taildeg_everytrial(i,a(find(aa==1)))=0;%去干扰
    % end
    a=taildeg_everytrial(i,:);
    for j=241:length(a)-240
        b=taildeg_everytrial(i,:);
        if (b(j)==0 && b(j+1)<0 && b(j-1)==0)||(b(j)==0 && b(j+1)>0 && b(j-1)==0)
            a(j:j+bin1)=0;
        elseif (b(j+1)==0 && b(j)<0 && b(j+2)==0)||(b(j+1)==0 && b(j)>0 && b(j+2)==0)
            a(j-bin2:j)=0;
        end
    end
    taildeg_everytrial1(i,:)=a;
end