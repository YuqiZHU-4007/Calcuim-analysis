%
clc;clear all;
[fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=10;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num)
T_spon={};T_US={};l_spon={};l_US={};T_non_spon=[];
for ii=1:trial.acq_block_num+1
    switch ii
        case 1
            t=(trial.spon_bef(2)-1)*frameb.per_cycle+1:(trial.spon_bef(3))*frameb.per_cycle;
            l_spon{ii}='Bef.';
        case mat2cell([2:trial.acq_block_num],1,ones(1,trial.acq_block_num-1))
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)]; 
        case trial.acq_block_num+1
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
            %l_spon{ii}='Test';
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)]; 
        case trial.acq_block_num+2
            t=(trial.test(3))*frameb.per_cycle+1:(trial.total)*frameb.per_cycle;
            l_spon{ii}='Aft.';
    end
    T_spon{ii}=t;
end
% us
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frameb.per_cycle+1:(trial.hab(3))*frameb.per_cycle;
            l_US{ii}='Hab';
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle ;
            
            l_US{ii}=['Acq.' num2str(ii-1)];
        case trial.acq_block_num+2
            t=(trial.test(2)-1)*frameb.per_cycle+1:(trial.test(3))*frameb.per_cycle;
            l_US{ii}='Tst';
    end
     T_US{ii}=t;
     T_non_spon=[T_non_spon,t];
end
stimUS=zeros(1,frameb.per_cycle*trial.total);ind_US=[];kk=1;
for ss=1:trial.acq_block_num
    for tt=1:trial.acq_block_trial
        ind=(trial.hab(3)+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frameb.per_cycle+frameb.us_start;
        stimUS(ind:ind+max(round(frameb.us_dur(ss)),1)-1)=1;%23
        ind_US(kk,1)=ind;kk=kk+1;
    end
end
trial_ind_CS=[trial.hab(2):trial.hab(3)];
for ii=1:trial.acq_block_num
    trial_ind_CS=[trial_ind_CS,[trial.hab(3)+1:trial.hab(3)+trial.acq_block_trial]+(ii-1)*(trial.acq_block_interval+trial.acq_block_trial)];
end
trial_ind_CS=[trial_ind_CS,trial.test(2):trial.test(3)];
stimCS=zeros(1,frameb.per_cycle*trial.total);ind_CS=[];kk=1;
for tt=trial_ind_CS
    stimCS((tt-1)*frameb.per_cycle+frameb.cs_start:(tt-1)*frameb.per_cycle+frameb.cs_end)=1;%23
    ind_CS(kk,1)=(tt-1)*frameb.per_cycle+frameb.cs_start;
    ind_CS(kk,2)=(tt-1)*frameb.per_cycle+frameb.cs_end;kk=kk+1;
end
ind_CS_hab=ind_CS(1:6,:);ind_CS_tst=ind_CS(end-5:end,:);
%load
path=uigetdir('H:\1.Test US\2.Tail free¡ª¡ªData from 117');
seg={'20201019','20201020','20201021','20201025','20201026','20201027'};
n=[3,3,4,4,4,3];kk=1;
Path=string;Name=string;
for ii=1:length(n)
    for jj=1:n(ii)
        name=fullfile(path,seg{ii},['fish' num2str(jj)],'Results_of_alltheta.mat');
        Path(kk,1)=name;
        Name(kk,1)=fullfile(seg{ii},['fish' num2str(jj)]);
        kk=kk+1;
    end
end
load(Path(5));name=fieldnames(Results_of_alltheta);
Data=struct;CR_trial={};V_trial_CS={};CR_loc={};V_trial_base={};V_trial_increased={};
for ii=1:length(name)
    Data=setfield(Data,name{ii},[]);
end
for zz=1:length(Path)
    if exist(Path(zz))
        disp(Path(zz))
        load(Path(zz));ind=find(strcmp(Path,Path(zz)));
        for ii=1:length(name)
            if min(size(getfield(Results_of_alltheta,name{ii})))<=1
                b=getfield(Data,name{ii})';a=getfield(Results_of_alltheta,name{ii});if ~isempty(a) b(ind,1:length(a))=a;end
            elseif min(size(getfield(Results_of_alltheta,name{ii})))>1
                b=getfield(Data,name{ii})';
                for jj=1:trial.acq_block_trial
                    bb=getfield(Results_of_alltheta,name{ii});
                    b{ind,jj}=bb(:,jj,:);
                end
            elseif min(size(getfield(Results_of_alltheta,name{ii})))>1 && strcmp(name{ii},'spon_env_loc')
                b=getfield(Data,name{ii})';
                for jj=1:size(getfield(Results_of_alltheta,name{ii}),1)
                    bb=getfield(Results_of_alltheta,name{ii});
                    b{ind,jj}=bb(jj,:,:);
                end
            end
            Data=setfield(Data,name{ii},b');
        end
        for hh=3
            switch hh
                case 1
                    ind_cs=ind_CS_hab;
                case 2
                    ind_cs=ind_CS_tst;
                case 3
                    ind_cs=ind_CS;
            end
            aa=(repmat(T_non_spon(env.env_loc(:,1)),length(ind_cs(:,1)),1)-repmat(ind_cs(:,1),1,length(env.env_loc(:,1))));
            ind=find(aa>=-0.1*60 & abs(aa)<=(frameb.cs_end-frameb.cs_start)+1-0.2*60);[~,ind]=ind2sub([size(aa,1),size(aa,2)],ind);
            env_loc_cs=[env.env_loc(ind,:),env.env_end(ind,:)];env_num=length(ind);
            %             [C]=mod(env.env_loc(:,1),frameb.per_cycle);ind=find(C<frameb.cs_start & C>=2*frameb.cs_start-(frameb.cs_end));
            %             env_loc_pre_cs=[env.env_loc(ind,:),env.env_end(ind,:)];
            %             figure,plot(alltheta(T_non_spon));hold on;plot(stimCS(T_non_spon)*180);
            %             scatter(env_loc_cs(:,1),alltheta(T_non_spon(env_loc_cs(:,1))),'r','filled');hold on
            %             scatter(env_loc_pre_cs(:,1),alltheta(T_non_spon(env_loc_pre_cs(:,1))),'k','filled');hold on     
            alltheta_n=normalize(alltheta,'zscore');%alltheta_n=normalize(alltheta,'range');
            dif=diff([alltheta_n(1);alltheta_n]);V=[];
            for ii=1:length(ind)
                V= [V;sum(abs(dif(T_non_spon(env_loc_cs(ii,1)):T_non_spon(env_loc_cs(ii,2)))))];
            end
            env_loc_pre_cs=[];V_base=[];
            for ii=1:length(T_spon)
                [~,env_loc,env_end,~]=env_detect(alltheta(T_spon{ii}),3*60,1);
                %h=figure;plot(alltheta(T_spon{ii}));hold on;scatter(env.env_loc(:,1),alltheta(T_spon{ii}(env.env_loc(:,1))),'filled');
                ind=T_spon{ii}([env_loc,env_end]);
                env_loc_pre_cs=[env_loc_pre_cs;ind];V_base_a=[];
                if ~isempty(ind)
                    for jj=1:size(ind,1)
                        V_base_a= [V_base_a;sum(abs(dif(ind(jj,1):ind(jj,2))))];
                    end
                else
                    V_base_a=[0 0];
                end
                V_base(:,ii)=sum(V_base_a)/length(T_spon{ii})*(frameb.cs_end-frameb.cs_start);
            end
            V_base=kron(V_base(1:trial.acq_block_num+1)',ones(trial.acq_block_trial,1));
            V_base=[repmat(V_base(1,:),trial.hab(1),1);V_base;repmat(V_base(end,:),trial.hab(1)-trial.acq_block_trial,1)];
            V_trial_base{hh}(:,zz)=V_base;
            [C,ia,~]=unique(fix(env_loc_cs(:,1)/frameb.per_cycle)+1);
            for ii=C'
                ind=find(fix(env_loc_cs(:,1)/frameb.per_cycle)+1==ii);
                V_trial_CS{hh}(ii,zz)=sum(V(ind));
                CR_trial{hh}(ii,zz)=length(ind);
                CR_loc{hh}(ii,zz)=mod(T_non_spon(env_loc_cs(ind(1),1)),frameb.per_cycle);
            end
        end
    end
end
V_trial_increased=V_trial_CS{3}./V_trial_base{3};
V_trial_increased_hab=mean(V_trial_increased(1:trial.hab(1),:),1);


learned=[2];
notsure=[];
acqlearned=[8 10];
nonlearned=[4 5];
all=[6,7,5,16,3,14,4,1,12,2,13];%[learned,notsure,acqlearned,nonlearned];
type=all;n=Name(type);

%CR amp.
b= CR_trial{3}(:,type);%b(find(b>=1))=1;n=Name(type);

a=[];
a(find(a>=1))=1;
aa=[];
aa(1,:)=mean(a(1:6,:));
aa(2,:)=mean(a(end-6+1:end,:));
[h,p]=ttest(aa(1,:),aa(2,:),'tail','left')