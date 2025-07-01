%%%%%mov-detection
function [taildeg_everytrial]=mov_detection_2(taildeg_everytrial,win,base1,base2,trailnum)
l=length(win);
for i=1:1:trailnum
    bin=win;
    a=taildeg_everytrial(i,:);
    for j=l:length(a)
        aa=sum(abs(diff(a(bin))))/l;
        if (sum(a(bin)<=base1)==l && a(j)<base2==1)%(aa<=base1 && a(j)<base2)
            taildeg_everytrial(i,j)=0;
        end
        bin=bin+1;
        if bin(end)>length(a)
            break;
        end
    end
    %taildeg_everytrial(i,find(diff(diff(taildeg_everytrial(i,:)))==0))=0;
end
