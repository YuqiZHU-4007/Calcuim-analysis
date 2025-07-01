function behavior=getmulti_motor_regresspr(re_startpoint_sd,startpoint_sd,y_3sd,frameb,fs)
%load(behavpath);
colorCS=[0.5 0.5 0.8];
if ~exist('re_start_end_point')
    re_start_end_point=[];
    re_start_end_point(:,1:2)=re_startpoint_sd(:,1:2);
    for ii=1:size(re_startpoint_sd,1)-1
        end_ind=find(y_3sd(startpoint_sd(ii)+1:startpoint_sd(ii+1)-1)~=0);
        end_ind=end_ind(end)+startpoint_sd(ii);
        re_start_end_point(ii,3)=end_ind;mod(end_ind,frameb.per_cycle);
    end
    end_ind=find(y_3sd(startpoint_sd(ii)+1:end)~=0);
    end_ind=end_ind(end)+startpoint_sd(ii);
    re_start_end_point(end,3)=end_ind;mod(end_ind,frameb.per_cycle);
%     end_ind=find(abs(y_3sd-[y_3sd(1),y_3sd(1:end-1)])~=0 & abs([y_3sd(2:end),0]-y_3sd)==0)';
%     end_ind=mod(end_ind,frameb.per_cycle);
%     re_start_end_point(:,3)=re_startpoint_sd(:,2)+32;
    %figure,plot(y_3sd);hold on;scatter(end_ind,y_3sd(end_ind));hold on;scatter(startpoint_sd,y_3sd(startpoint_sd),'filled');
end
% motor_CS=re_start_end_point(find(re_start_end_point(:,2,:)>= frameb.cs_start & re_start_end_point(:,2,:)< frameb.cs_end),:);
% motor_spon=re_start_end_point(find(re_start_end_point(:,2,:)< frameb.cs_start | re_start_end_point(:,2,:)>= frameb.cs_end),:);
motor_CS=(find(re_start_end_point(:,2,:)>= frameb.cs_start & re_start_end_point(:,2,:)< frameb.cs_end));
motor_spon=(find(re_start_end_point(:,2,:)< frameb.cs_start | re_start_end_point(:,2,:)>= frameb.cs_end));
time_be=0:1/fs.behavior:(length(y_3sd)-1)/fs.behavior;
behavior=zeros(5,length(time_be));
for hh=1:3
    switch hh
        case 1 %全长，所有shaking
            y=y_3sd; %figure,plot(abs(y));
        case 2 %only CS
            y= zeros(size(y_3sd));
            ind=[startpoint_sd(motor_CS);re_start_end_point(motor_CS,3)']';%[(motor_CS(:,1)-1)*frameb.per_cycle+motor_CS(:,2) (motor_CS(:,1)-1)*frameb.per_cycle+motor_CS(:,3)];
            for jj=1:size(ind,1)
                y(ind(jj,1):ind(jj,2))=y_3sd(ind(jj,1):ind(jj,2));
            end
            %             figure,plot(abs(y));hold on;
            %             a=[0 20];
            %             patch1=patch([[frameb.cs_start:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_start:frameb.per_cycle:size(y',1)]']',...
            %                 repmat([min(a) min(a) max(a) max(a)],length([frameb.cs_start:frameb.per_cycle:size(y',1)]),1)',...
            %                 colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        case 3 %only Spon
            y= zeros(size(y_3sd));
            ind=[startpoint_sd(motor_spon);re_start_end_point(motor_spon,3)']';%[(motor_spon(:,1)-1)*frameb.per_cycle+motor_spon(:,2) (motor_spon(:,1)-1)*frameb.per_cycle+motor_spon(:,3)];
            for jj=1:size(ind,1)
                y(ind(jj,1):ind(jj,2))=y_3sd(ind(jj,1):ind(jj,2));
            end
            %             figure,plot(abs(y));hold on;
            %             a=[0 20];
            %             patch1=patch([[frameb.cs_start:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
            %                 [frameb.cs_start:frameb.per_cycle:size(y',1)]']',...
            %                 repmat([min(a) min(a) max(a) max(a)],length([frameb.cs_start:frameb.per_cycle:size(y',1)]),1)',...
            %                 colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    end
    %     time_ca=0:fs.ca:(length(activities_preCS_dfdf_aftcorrect(:,1))-1)*fs.ca;
    behavior(hh,:)=abs(y);
end
figure,plot(behavior(1,:),'k');hold on;
plot(behavior(2,:),'r');hold on;
plot(behavior(3,:),'b');hold on;
legend('all','CS','Spon.');
a=[0 20];
patch1=patch([[frameb.cs_start:frameb.per_cycle:size(y',1)]'...
    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
    [frameb.cs_end:frameb.per_cycle:size(y',1)]'...
    [frameb.cs_start:frameb.per_cycle:size(y',1)]']',...
    repmat([min(a) min(a) max(a) max(a)],length([frameb.cs_start:frameb.per_cycle:size(y',1)]),1)',...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
%scatter(startpoint_sd+32,abs(y_3sd(startpoint_sd+32)));hold on;scatter(startpoint_sd,abs(y_3sd(startpoint_sd)),'filled');