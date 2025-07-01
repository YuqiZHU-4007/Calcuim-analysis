clc;
activities_new={};
[inputpath]=uigetdir('G:\data','activities.xlsx');
listdir = dir([inputpath '/*.xlsx']);
outputpath=uigetdir('G:\data','outputpath');
if ~isempty(listdir)
    for ii = 1:length(listdir)
        filename = listdir(ii).name;
        activities=xlsread([inputpath '\' filename]);
        activities_new{ii,1}=activities(:,3:4:size(activities,2));
    end
end
%不同的文件
% for i=1:n
% activities=xlsread([inputpath 'Z' num2str(i,'%02d') '\'  'Z' num2str(i,'%02d')]);
% activities_new{i,1}=activities(:,3:4:size(activities,2));
% end

% %不同sheet
% [~,sheet,~]=xlsfinfo([inputpath,inputname]);
% for i=1:size(sheet,2)
% activities=xlsread([inputpath,inputname],sheet{i});
% activities_new{i,1}=activities(:,3:4:size(activities,2));
% end
% 
% 
save([outputpath '\activities_new.mat'],'activities_new');