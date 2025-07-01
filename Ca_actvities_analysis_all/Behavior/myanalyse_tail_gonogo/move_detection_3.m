function [t_dur_cell,taildeg_everytrial]=move_detection_3(t_dur,taildeg_everytrial,trailnum,base1,base2,base3,bin2)
for i=1:trailnum
    a=t_dur(:,1,i);%a(find(a==0),:)=[];
    b=t_dur(:,2,i);%b(find(b==0),:)=[];
    for j=1:length(a)
        ind=zeros;
        if a(j)~=0
            text=taildeg_everytrial(i,a(j):a(j)+b(j)-1).*10;
            dif=abs(diff(text));
            ind=find(abs(diff(text))>=50);%%50 认为是高峰 与低尾动干扰的阈值
            if (sum(abs(text))/length(text)<=base1 && length(ind)==0 ) %%均值小于阈值，且持续时间长于120&& b(j)>=120 判断正负干扰全在一侧
                t_dur(j,:,i)=0;taildeg_everytrial(i,a(j):a(j)+b(j)-1)=0;%%去干扰
            elseif (dif(1)<=base2  && length(ind)~=0) %%%去峰值前端
                ind=find([dif(2:end) dif(end)]>=base2 &  dif>=base2 & [dif(1) dif(1:end-1)]<base2);%[dif(3:end) dif(end) dif(end)]>=base2 &
                %taildeg_everytrial1(i,a(j):a(j)+ind(1))=0;
                if length(ind)~=0
                    ind(1)=ind(1)+a(j);
                    tind=detect_tail(taildeg_everytrial(i,:).*10,ind(1),a(j)+b(j)-1,base2,base3);
                    t_dur(j,:,i)=[ind(1) tind];
                    taildeg_everytrial(i,a(j):ind(1)-1)=0;
                end
            end
            if (b(j)>200 && length(ind)~=0) %%%去双峰之间
%                 if (sum(abs(text(50:80)))/length(text(50:80))<=base1 && length(find(abs(diff(text(50:80)))>=50))==0)%%如果有很长的连续动，加这个限定条件
                    ind=find(abs(diff(text))>=45);
                    mind=find(ind-[ind(2:end) ind(end)]<=-base3);%%%峰值连续变化index的插值的距离，大于30 认为有第二次动
                    %mind1=find([ind(2:end) ind(end)]-ind>=base3);
                    l=length(mind);
                    if l~=0
                        mind1=ind(mind)+a(j);mind=ind(mind+1)+a(j);
                        tind=detect_tail(taildeg_everytrial(i,:).*10,[t_dur(j,1,i) mind],[mind1 t_dur(j,1,i)+t_dur(j,2,i)-1],base2,base3);
                        %taildeg_everytrial(i,a(j)+tind(1)-1:mind(1)-1)=0;
                        t_dur(j,2,i)=tind(1);
                        t_dur(size(t_dur,1)+1:size(t_dur,1)+l,:,i)=[mind' tind(2:end)'];
                    end
                end
%             end
        end
    end
    tt=sortrows(t_dur(:,:,i),1);
    tt(find(tt(:,1)==0),:)=[];tt(:,3)=tt(:,1)+tt(:,2)-1;t_dur_cell{i,1}=tt;
end