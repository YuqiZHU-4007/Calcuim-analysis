%[inputpath]=uigetdir('E:\A_Data_lightsheet\Data_vmat\');
actpathes=batch_code('txt of <activties_new> path') ;
folder1='total';
for batchi=4:length(actpathes)
    actpath= actpathes{batchi};
    actpath=[actpath(1:46) 'subregoin_fromROI']
    listdir = dir([actpath '/*.*']);
    fishname=actpath(32:45)
    for ii=1:length(listdir)
        if ~strcmp(listdir(ii).name,'.') || ~strcmp(listdir(ii).name,'..')
            sublist=dir([listdir(ii).folder '\' listdir(ii).name '\' folder1]);
            str=[listdir(ii).name '-mean-raw trace.fig'];
            if sum(strcmp(str,{sublist.name}))==1
                state='on';listdir(ii).name
                pathind=find(strcmp(str,{sublist.name}));
                path=[sublist(pathind).folder '\' sublist(pathind).name];
                h=openfig(path);
                h.Visible='on';
                h.Name=listdir(ii).name;
                for zz=[4 1 1]
                %set(h. Children(zz),'visible','off');
                delete(h. Children(zz))
                end  
                set(h,'OuterPosition',[0 0 1 1],'position',[0.2,0.2,0.75,0.4]);
                set(h. Children(1),'position',[0.05,0.15,0.9,0.7]);
                h. Children(1).Title.String=[listdir(ii).name '-Acq.'];
%                 set(h,'OuterPosition',[0 0 1 1],'position',[0.2,0.2,0.4,0.6]);
%                 set(h. Children(3),'Position',[0.13 0.167763157894737 0.75551724137931 0.757236842105262]);
%                 %set(h. Children(3).Title.String,fishname);
%                 h. Children(3).Title.String=fishname;
%                 %set(h. Children(1),'visible','off')
%                 delete(h. Children(1))
%                 print('-clipboard','-dmeta')
                pause(state);
                a=input('state\n'); if a==1 state='off'; end
                pause(state);
                close(h);
            else
                continue;
            end
        else
            continue;
        end
    end
end

for batchi=4:length(actpathes)
    actpath= actpathes{batchi};
    actpath=[actpath(1:46) 'subregoin_fromROI']
    listdir = dir([actpath '/*.*']);
    fishname=actpath(32:45)
    for ii=1:length(listdir)
        if ~strcmp(listdir(ii).name,'.') || ~strcmp(listdir(ii).name,'..')
            sublist=dir([listdir(ii).folder '\' listdir(ii).name]);
            if sum(strcmp('trace_hab_test.fig',{sublist.name}))==1
                state='on';listdir(ii).name
                pathind=find(strcmp('trace_hab_test.fig',{sublist.name}));
                path=[listdir(ii).folder '\' listdir(ii).name '\' sublist(pathind).name];
                h=openfig(path);
                h.Visible='on';
                h.Name=listdir(ii).name;
                for zz=[1  1]
                    %set(h. Children(zz),'visible','off');
                    delete(h. Children(zz))
                end
                set(h,'OuterPosition',[0 0 1 1],'position',[0.2,0.2,0.75,0.4]);
                 set(h. Children(2),'position',[0.1,0.2,0.4,0.7]);
                 set(h. Children(1),'position',[0.6,0.2,0.2,0.7]);
                delete(h. Children(1).Children(1));
                delete(h. Children(1).Children(1));
                h. Children(2).Title.String=[listdir(ii).name '-Hab. and test'];
%                 set(h,'OuterPosition',[0 0 1 1],'position',[0.2,0.2,0.2,0.5]);
%                 set(h. Children(3),'Position',[0.218453188602442 0.11 0.680429790120962 0.815]);
%                 %set(h. Children(3).Title.String,fishname);
%                 h. Children(3).Title.String=fishname;
%                 C=h. Children;
%                 %set(h. Children(1),'visible','off')
%                 delete(C(1));delete(C(2));delete(C(4))
%                 print('-clipboard','-dmeta')
                pause(state);
                a=input('state\n'); if a==1 state='off'; end
                pause(state);
                close(h);
            else
                continue;
            end
        else
            continue;
        end
    end
end

para=batch_code('txt of <para> path') ;
beh=batch_code('txt of <behav> path') ;
max_amp=[];s_point={};
for batchi=1:length(para)
    parapath= para{batchi};
    behaviorpath=beh{batchi};
    load(parapath);
    a=load(behaviorpath);
    fishname=parapath(32:45)
%     [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(a.delta_r_bef1,re_startpoint,fishname,[],[]);
%     figure,plot_behavior_onset(a.delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
    startt=re_startpoint(find(re_startpoint(:,1)>=trial.test(2) & re_startpoint(:,1)<=trial.test(2)+3-1 & ...
        re_startpoint(:,2)<=frameb.cs_end & re_startpoint(:,2)>=frameb.cs_start),:);
    start=zeros(3,2);
    for zz=trial.test(2):trial.test(2)+3-1
        ind=find(startt(:,1)==zz);
        if length(ind)==1
            start(zz-trial.test(2)+1,:)=startt(ind,:);
        elseif length(ind)>1
            ind=ind(1);
        start(zz-trial.test(2)+1,:)=startt(ind,:);
        elseif isempty(ind)
            start(zz-trial.test(2)+1,:)=nan(1,2);
        end
    end
    for zz=1:size(start,1)
        Y=y_3sd(startpoint(find(re_startpoint(:,1)==start(zz,1) & re_startpoint(:,2)==start(zz,2))):end);
        if isempty(Y)
            max_amp(zz,batchi)=nan;
        else
        ind=find(Y==0 & [Y(2:end) Y(end)]==0);
        ind=ind(1);
        max_amp(zz,batchi)=max(abs(Y(1:ind)));
        end
    end
    start(:,2)=(start(:,2)-frameb.cs_start)/fs.behavior;
    s_point{batchi}=start;
    %a=input('state\n'); 
end
