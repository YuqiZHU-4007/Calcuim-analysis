clc;clear all;close all;
mat_name=[];
Path{1}={['H:\1.Test US\2.Tail free！！Data from 117\20220709\fish1\',mat_name],...%no Cb after Learning %['H:\1.Test US\2.Tail free！！Data from 117\20220802\fish2\',mat_name],...%not typical
    ['H:\1.Test US\2.Tail free！！Data from 117\20220803\fish2\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220814\fish1\',mat_name],...%['H:\1.Test US\2.Tail free！！Data from 117\20220915\fish1\',mat_name],...%calcium wrong
    ['H:\1.Test US\2.Tail free！！Data from 117\20210402\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210709\fish2\',mat_name]
    %     ['H:\1.Test US\2.Tail free！！Data from 117\20210805\fish1\',mat_name],...
    %     ['H:\1.Test US\2.Tail free！！Data from 117\20210430\fish1\',mat_name],...%SPACE 8*3 SKIP
    %     ['E:\A_Data_lightsheet\Data_huc\20190604\fish4\',mat_name],... %SPACE 5*6
    %     ['E:\A_Data_lightsheet\Data_huc\20190605\fish3\',mat_name],...%SPACE 5*6
    %     ['E:\A_Data_lightsheet\Data_huc\20190609\fish3\',mat_name],...%SPACE 5*6
    %     ['E:\A_Data_lightsheet\Data_huc\20190514\fish3\',mat_name] % MASSIVE
    };%learner
Path{2}={
    ['H:\1.Test US\2.Tail free！！Data from 117\20210514\fish2\',mat_name],...%un shuffle
    ['H:\1.Test US\2.Tail free！！Data from 117\20210702\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220831\fish1\',mat_name ]
    %['H:\1.Test US\2.Tail free！！Data from 117\20210802\fish2\',mat_name],...
    %     ['E:\A_Data_lightsheet\Data_huc\20190612\fish1\',mat_name ],...
    %     ['E:\A_Data_lightsheet\Data_huc\20190612\fish2\',mat_name ]
    %['H:\1.Test US\2.Tail free！！Data from 117\20210617\fish2\',mat_name],...
    };%unpair
Path{3}={['H:\1.Test US\2.Tail free！！Data from 117\20220803\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220816\fish2\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220824\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20220829\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210510\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210529\fish1\',mat_name]
    %     ['G:\data_huc_non_learner\20190604\fish2\',mat_name ],...
    %     ['G:\data_huc_non_learner\20190607\fish3\',mat_name ]
    };%non-learner
Path{4}={['H:\1.Test US\2.Tail free！！Data from 117\20220817\fish1\',mat_name ],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210401\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210717\fish1\',mat_name],...
    ['H:\1.Test US\2.Tail free！！Data from 117\20210822\fish1\',mat_name]
 %       ['H:\1.Test US\2.Tail free！！Data from 117\20210719\fish1\',mat_name],...
    };%faded-learner
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
save([savepath '\Path'],'Path','-v7.3');

pp=string;
for ii=1:length(Path)
    path=Path{ii};
    for jj=1:length(path)
        checkpath(path{jj});
        a=path{jj}(end-14:end);
        pp=cat(1,pp,a);
    end
end


labels_typei={'CS-avtivation','CS-inhibition','US-activation','CS-up regulate','CS-down regulate','CS-stable regulate','US-up regulate','US-down regulate','US-stable regulate'};
for ii=1:4
    path=Path{ii};
    for jj=1:length(path)
        file=dir([path{jj},'/*.jpg']);
        for zz=1:length(file)
        delete(fullfile(file(zz).folder,file(zz).name))
        end
        for zz=1:length(labels_typei)
         try 
             rmdir(string(fullfile(path{jj},labels_typei(zz))),'s')
         catch 
             warning(string(fullfile(path{jj},labels_typei(zz))))
         end
        end
    end
end

for ii=1:2
    path=Path{ii};
    for jj=1:length(path)
        US_mapping(path{jj})
    end
end
% Path{1}={%['H:\1.Test US\2.Tail free！！Data from 117\20210402\fish1\',mat_name ],
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210430\fish1\',mat_name],...
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210709\fish2\',mat_name],...
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210805\fish1\',mat_name],...
%     ['E:\A_Data_lightsheet\Data_huc\20190514\fish3\',mat_name],...
%     ['E:\A_Data_lightsheet\Data_huc\20190604\fish4\',mat_name],...
%     ['E:\A_Data_lightsheet\Data_huc\20190609\fish3\',mat_name]
%     };%learner
% Path{2}={%['H:\1.Test US\2.Tail free！！Data from 117\20210514\fish2\',mat_name ],...
%     %['H:\1.Test US\2.Tail free！！Data from 117\20210617\fish2\',mat_name],...
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210702\fish1\',mat_name],...
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210802\fish2\',mat_name],...
%     ['E:\A_Data_lightsheet\Data_huc\20190612\fish1\',mat_name ],...
%     ['E:\A_Data_lightsheet\Data_huc\20190612\fish2\',mat_name ]
%     };%unpair
% Path{3}={['H:\1.Test US\2.Tail free！！Data from 117\20210510\fish1\',mat_name ],...
%     ['H:\1.Test US\2.Tail free！！Data from 117\20210529\fish1\',mat_name],...
%     ['G:\data_huc_non_learner\20190604\fish2\',mat_name ],...
%     ['G:\data_huc_non_learner\20190607\fish3\',mat_name ]
%     };%non-learner
% %Path{4}={['H:\1.Test US\2.Tail free！！Data from 117\20210719\fish1\',mat_name ],...
%     %['H:\1.Test US\2.Tail free！！Data from 117\20210717\fish1\',mat_name]};%fadd-learner

