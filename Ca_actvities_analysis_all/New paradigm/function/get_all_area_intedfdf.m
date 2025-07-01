%[inputpath]=uigetdir('E:\A_Data_lightsheet\Data_vmat\');
actpathes_dim=batch_code('txt of <activties_new_dim> path') ;
locpathes_dim=batch_code('txt of <location_dim> path') ;
parapathes_dim=batch_code('txt of <para_dim> path') ;

% actpathes_fla=batch_code('txt of <activties_new_fla> path') ;
% locpathes_fla=batch_code('txt of <location_fla> path') ;
% parapathes_fla=batch_code('txt of <para_fla> path') ;

% actpathes=[actpathes_dim;actpathes_fla];
% locpathes=[locpathes_dim;locpathes_fla];
% parapathes=[parapathes_dim;parapathes_fla];
actpathes=[actpathes_dim];
locpathes=[locpathes_dim];
parapathes=[parapathes_dim];

actpathes={actpathes_dim{2:end}};
locpathes={locpathes_dim{2:end}};
parapathes={parapathes_dim{2:end}};

file_str='activities_aft_process.mat';
zz=0;CS_acq_block=[];CS_acq_ind=[];CS_hab_tst=[];CS_hab_tst5=[];CS_acq_acti=[]; CS_acti=[];
for batchi=1:length(actpathes)
    actpath= actpathes{batchi};
    locpath= locpathes{batchi};
    parapath= parapathes{batchi};
    load(parapath)
    listdir = dir([actpath(1:end-18) '/*.*']);
    fishname=actpath(32:45)
    for jj=1:length(listdir)
        if ~strcmp(listdir(jj).name,'.') || ~strcmp(listdir(jj).name,'..')
            if sum(strcmp(file_str,{listdir(jj).name}))==1
                state='on';listdir(jj).name
                pathind=find(strcmp(file_str,{listdir(jj).name}));
                path=[listdir(jj).folder '\' listdir(jj).name];
                a=load(path);
                cut_move=true;
                ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
                [rawacti,area]=calculate_integtate_dfdf_main(a.activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
                [~,area_hab5]=calculate_integtate_dfdf_main_hab5(a.activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
                zz=zz+1;
                %extra raw_acti and integrated dfdf
                [~,sheet,~]=xlsfinfo(locpath);
                num=[];
                for kk=1:size(sheet,2)
                    num=xlsread(locpath,sheet{kk});
                    if sum(strcmp('Pr-5HT',{sheet{kk}}))==1
                        sheet{kk}='Pr5HT';
                    elseif sum(strcmp('Pr-DA',{sheet{kk}}))==1
                        sheet{kk}='PrDA';
                    elseif sum(strcmp('Sub pallium',{sheet{kk}}))==1
                        sheet{kk}='SP';
                        elseif sum(strcmp('CG',{sheet{kk}}))==1
                        sheet{kk}='GC';
                    end
                    if ~isfield(CS_acq_block,sheet{kk})
                        CS_acq_block=setfield(CS_acq_block,sheet{kk},[]);
                        CS_hab_tst=setfield(CS_hab_tst,sheet{kk},[]);
                        CS_hab_tst5=setfield(CS_hab_tst5,sheet{kk},[]);
                        CS_acq_ind=setfield(CS_acq_ind,sheet{kk},[]);
                        CS_acq_acti=setfield(CS_acq_acti,sheet{kk},[]);
                         CS_acti=setfield( CS_acti,sheet{kk},[]);
                    end
                    disp(sheet{kk});
                    area_num_acq=[];ind=[];area_num_hab=[]; area_num_hab5=[];acti_acqa=[];acti_acq=[];acti_all=[];
                    for ii=1:size(num,1)
                        acti=a.activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
                        acti_acqa=[];
                        for jjj=trial.acq(2):trial.acq(3)
                            acti_acqa=[acti_acqa;acti((jjj-1)*frame.per_cycle+frame.cs_start:(jjj-1)*frame.per_cycle+frame.us_start-1,:)];
                        end
                        %rawacti_acq_mean=rawacti.acq_mean_blocks{num(ii,1),1}{num(ii,1),1};
                        acti_acq(ii,1:length(acti_acqa))=acti_acqa';
                        acti_all(ii,1:length(acti(trial.hab(3)*frame.per_cycle+1:trial.acq(3)*frame.per_cycle)))=acti(trial.hab(3)*frame.per_cycle+1:trial.acq(3)*frame.per_cycle);
                        area_num_acq(ii,1:length(area.CS_acq_block{num(ii,1),1}(:,num(ii,2))'))=area.CS_acq_block{num(ii,1),1}(:,num(ii,2))';
                        area_num_hab(ii,1:length(area.CS_hab_tst{num(ii,1),1}(:,num(ii,2))'))=area.CS_hab_tst{num(ii,1),1}(:,num(ii,2))';
                        area_num_hab5(ii,1:length(area_hab5.CS_hab_tst{num(ii,1),1}(:,num(ii,2))'))=area_hab5.CS_hab_tst{num(ii,1),1}(:,num(ii,2))';
                        %              area_num_acq_1st(:,ii)=a.area.CS_acq_1st_block{num(ii,1),1}(:,num(ii,2));
                        %              area_num_hab_tst(:,ii)=a.area.CS_hab_tst{num(ii,1),1}(:,num(ii,2));
                        %a_event_rep=[a_event_rep a_event{ii,1}.ind];
                        ind(ii,:)=zz;
                    end
                    CS_acq_block=setfield( CS_acq_block,sheet{kk},...
                        [[getfield(CS_acq_block,sheet{kk}),zeros(size(getfield(CS_acq_block,sheet{kk}),1),size(area_num_acq,2)-size(getfield(CS_acq_block,sheet{kk}),2))];...
                        area_num_acq]);
                    CS_hab_tst=setfield( CS_hab_tst,sheet{kk},...
                        [[getfield(CS_hab_tst,sheet{kk}),zeros(size(getfield(CS_hab_tst,sheet{kk}),1),size(area_num_hab,2)-size(getfield(CS_hab_tst,sheet{kk}),2))];...
                        area_num_hab]);
                    CS_hab_tst5=setfield( CS_hab_tst5,sheet{kk},...
                        [[getfield(CS_hab_tst5,sheet{kk}),zeros(size(getfield(CS_hab_tst5,sheet{kk}),1),size(area_num_hab5,2)-size(getfield(CS_hab_tst5,sheet{kk}),2))];...
                        area_num_hab5]);
                    CS_acq_ind=setfield( CS_acq_ind,sheet{kk},...
                        [getfield(CS_acq_ind,sheet{kk});ind]);
                    CS_acq_acti=setfield(CS_acq_acti,sheet{kk},...
                        [getfield(CS_acq_acti,sheet{kk});acti_acq]);%,zeros(size(getfield(CS_acq_acti,sheet{kk}),1),size(acti_acq,2)-size(getfield(CS_acq_acti,sheet{kk}),2))];...  
                     CS_acti=setfield(CS_acti,sheet{kk},...
                        [getfield(CS_acti,sheet{kk});acti_all]);
                end
            end
        end
    end
end

name=fieldnames(CS_acq_block)
X1=[];X2=[];X3=[];XX1=[];XX2=[];XX3=[];X4=[];X5=[];
m1=[];sd1=[];m=[];sd=[];n=[];
for ii=1:length(name)
    areaa=getfield(CS_acq_block,name{ii})';
    ind=getfield(CS_acq_ind,name{ii})';ind=unique(ind);
    win=1:size(areaa,1);
    if ~isempty(areaa)
    m=mean(areaa(win,:),2);
    sd=std(areaa(win,:),[],2);
    n=size(areaa(win,:),2);
    X1=setfield(X1,name{ii},[m,sd,n*ones(size(m))]);
    XX1=setfield(XX1,name{ii},areaa(win,:));
    else
         X1=setfield(X1,name{ii},[]);
         XX1=setfield(XX1,name{ii},[]);
    end
    
    areaa=getfield(CS_hab_tst,name{ii})';
    if ~isempty(areaa)
    m1=mean(areaa,2);
    sd1=std(areaa,[],2);
    %X2=setfield(X2,name{ii},[[m1(1,ii);m(:,ii)],[sd1(1,ii);sd(:,ii)],n*ones(size([m1(1,ii);m(:,ii)]))]);
    X2=setfield(X2,name{ii},[m1,sd1,n*ones(size(m1))]);
    XX2=setfield(XX2,name{ii},areaa);
    else
         X2=setfield(X2,name{ii},[]);
         XX2=setfield(XX2,name{ii},[]);
    end
    
    areaa=getfield(CS_hab_tst5,name{ii})';
    if ~isempty(areaa)
    m1=mean(areaa,2);
    sd1=std(areaa,[],2);
    X3=setfield(X3,name{ii},[[m1(1,:);m],[sd1(1,:);sd],n*ones(size([m1(1,:);m]))]);
     XX3=setfield(XX3,name{ii},areaa);
     else
         X3=setfield(X3,name{ii},[]);
          XX3=setfield(XX3,name{ii},[]);
    end
    
    areaa=getfield(CS_acq_acti,name{ii})';
    win=1:size(areaa,1);
    if ~isempty(areaa)
    X4=setfield(X4,name{ii},areaa);
    else
         X4=setfield(X4,name{ii},[]);
    end
    
    areaa=getfield(CS_acti,name{ii})';
    win=1:size(areaa,1);
    if ~isempty(areaa)
    X5=setfield(X5,name{ii},areaa);
    else
         X5=setfield(X5,name{ii},[]);
    end
end

%%acq
X=[];m1=[];sd1=[];
for jj=1:length(name)
    for ii=1:zz
        ind=find(getfield(CS_acq_ind,name{jj})==ii);
        areaa=getfield(CS_acq_block,name{jj})';
        areaa=areaa(:,ind);
        m1(:,ii)=mean(areaa,2);
        sd1(:,ii)=std(areaa,[],2);
    end
    X=setfield(X,name{jj},m1);
end
X_acq1=[];X_acq2=[];
name_grap={'Raphe','LC','GC','PrDA','Pr5HT','PT','Pallium','SP','OB'};
for jj=1:length(name_grap)
    a=getfield(X,name_grap{jj});
    X_acq1(jj,:)=a(1,:);
    X_acq2(jj,:)=a(end,:);
end

%%hab_tst
X=[];m1=[];sd1=[];
for jj=1:length(name)
    for ii=1:zz
        ind=find(getfield(CS_acq_ind,name{jj})==ii);
        areaa=getfield(CS_hab_tst,name{jj})';
        areaa=areaa(:,ind);
        m1(:,ii)=mean(areaa,2);
        sd1(:,ii)=std(areaa,[],2);
    end
    X=setfield(X,name{jj},m1);
end
X_bef=[];X_aft=[];
name_grap={'Raphe','LC','GC','PrDA','Pr5HT','PT','Pallium','SP','OB'};
for jj=1:length(name_grap)
    a=getfield(X,name_grap{jj});
    X_bef(jj,:)=a(1,:);
    X_aft(jj,:)=a(end,:);
end
