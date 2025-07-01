function activities_new=repmat_act_to_cell_format(act_path,env_path)

load(act_path);
load(env_path);
T=opt.T;
T=T-1;
activities_new=cell(T,1);

for ii=1:T
    ind=find(env.supervoxel(:,3)==ii);
    activities_new{ii,1}=activities(ind,:)';
end


if isdir(act_path)
save([act_path '\activities_new.mat'],'activities_new');
else
     [filepath,~,~] = fileparts(act_path);
     save([filepath '\activities_new.mat'],'activities_new');
end
end