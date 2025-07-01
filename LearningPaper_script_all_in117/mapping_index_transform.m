clc;clear all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
group={'Learner','Unpair','Non-learner','Faded-learner'};
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
labels_type=load(fullfile(fullfile(Path{2}{1},'singletrial'),'CR_ind_summary_singletrial.mat'), 'labels_typei');
columns = {'t','x', 'y', 'z'};
align_label={'All','Go','NoGo'};
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'};
for ind_type=1:2
    switch ind_type
        case 1
            savepath=checkpath('W:\ind_all_CSUS_RESPONSIVE');
        case 2
            savepath=checkpath('W:\ind_all_emergeCSUS_inthissession');
        case 3
            savepath=checkpath('W:\ind_all_emergeCSUS_precond');
        case 4
            savepath=checkpath('W:\ind_all_emergeCSUS_postcond');
    end
    for batchi=1:4
        path=Path{batchi};align_win_total=zeros(length(sessionx),3);
        for typei=1:length(labels_type.labels_typei)
            sap=checkpath(fullfile(savepath,group{batchi},string(labels_type.labels_typei(typei))));
            for sessioni=1:length(sessionx)
                All=cell(1,3);
                for aligntypei=1:3
                    for fishi=num{batchi}
                        pp=path{fishi};p=checkpath(fullfile(pp,'singletrial'));nn=[pp(end-14:end-7),pp(end-5:end-1)];
                        load(fullfile(pp, '/activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo');
                        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                        load(fullfile(p,'CR_ind_summary_singletrial.mat'), 'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS_inthissession');
                        typeind=find(strcmp(labels_type.labels_typei(typei),labels_typei));
                        switch aligntypei
                            case 1
                                win=align_win;
                            case 2
                                win=align_win_go;
                            case 3
                                win=align_win_nogo;
                        end
                        ind=win(:,sessioni);ind(isnan(ind))=[];
                        align_win_total(sessioni,aligntypei)=align_win_total(sessioni,aligntypei)+length(ind);
                        for indi=1:length(ind)
                            switch ind_type
                                case 1
                                    a=ind_all_CSUS_RESPONSIVE{typeind,ind(indi),1};
                                case 2
                                    a=ind_all_emergeCSUS_inthissession{typeind,ind(indi),1,1};
                                case 3
                                    a=ind_all_emergeCSUS_inthissession{typeind,ind(indi),1,2};
                                case 4
                                    a=ind_all_emergeCSUS_inthissession{typeind,ind(indi),1,3};
                            end
                            %%ind_all_CSUS_RESPONSIVE{typeind,ind(indi),1};
                            All{1,aligntypei}=cat(1,All{1,aligntypei},supervolxeli(a,1:3));
                        end
                    end
                    if ~isempty(All{1,aligntypei})
                        t=ones(size(All{1,aligntypei},1),1);
                        M=table(t,All{1,aligntypei}(:,1),All{1,aligntypei}(:,2),All{1,aligntypei}(:,3),'VariableNames', columns);
                        name=strcat(labels_type.labels_typei(typei),'_',sessionx(sessioni),'_',string(align_label{aligntypei}),'.txt');
                        outpath=fullfile(sap,name);
                        writetable(M,outpath,'Delimiter','space','WriteVariableNames',false);
                    end
                end
            end
        end
        save([fullfile(savepath,group{batchi}),'\align_win_total.mat'],'align_win_total');
    end
end