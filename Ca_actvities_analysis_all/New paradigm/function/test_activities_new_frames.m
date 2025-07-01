%test activities_new «∑Ò…Ÿ÷°
actpathes=batch_code('txt of <activties_new> path') ;
para_pathes=batch_code('txt of <para> path') ;

nbatch = length(actpathes);

for batchi=1:nbatch
    actpath= actpathes{batchi};
    load(actpath);
    parapath=para_pathes{batchi};
    load(parapath);
    %input('1');
end

activities={};
for ii = 1:size(activities_new,1)
    if ~isempty(activities_new{ii,1}(:))
        activities{ii,1}=[activities_new{ii,1}(1,:);activities_new{ii,1};activities_new{ii,1}(size(activities_new{ii,1},1),:);];
    else
        activities{ii,1}=[];
    end
    %activities{ii,1}(trial.acq(3)*frame.per_cycle+1:(trial.acq(3)+2)*frame.per_cycle,:)=[];
end
for jj=trial.acq(2)+2*trial.acq_block_trial:trial.acq(2)+2*trial.acq_block_trial+3-1
    for ii = 1:size(activities_new,1)
        activities{ii,1}((jj-1)*frame.per_cycle+1:end+1,:)=[activities{ii,1}((jj-1)*frame.per_cycle+1,:);...
            activities{ii,1}((jj-1)*frame.per_cycle+1:end,:)];
    end
end
for ii = 1:size(activities_new,1)
    if ~isempty(activities_new{ii,1}(:))
        activities_new{ii,1}= activities{ii,1}(2:end-1,:);
    else
        activities_new{ii,1}=[];
    end
    %activities{ii,1}(trial.acq(3)*frame.per_cycle+1:(trial.acq(3)+2)*frame.per_cycle,:)=[];
end
figure,plot_rawtrace_trials(activities{1,1}(:,1),[],fs,frame,trial,startpoint,1);
save([actpath],'activities_new','-v7.3');
% for jj=trial.test(2)
%     for ii = 4:7
%         activities_new{ii,1}((jj-1)*frame.per_cycle+1:end+1,:)=[activities_new{ii,1}((jj-1)*frame.per_cycle+1,:);...
%             activities_new{ii,1}((jj-1)*frame.per_cycle+1:end,:)];
%     end
% end
% 
[behaviorname,behaviorpath]=uigetfile(['E:\A_Datlightsheet\Data_vmat\' '.mat'],'behavior');
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([behaviorpath behaviorname]);
save([behaviorpath '\para'],'fs','time','frame','frameb','trial','re_startpoint','startpoint','y_3sd','ishuc','-v7.3');