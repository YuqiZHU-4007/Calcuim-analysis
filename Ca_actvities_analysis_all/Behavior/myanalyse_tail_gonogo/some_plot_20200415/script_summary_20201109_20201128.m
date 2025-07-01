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
seg={'20201109','20201110','20201111','20201126','20201127','20201128'};
n=[4,5,6,4,5,3];kk=1;
Path=string;Name=string;
for ii=1:length(n)
    for jj=1:n(ii)
        name=fullfile(path,seg{ii},['fish' num2str(jj)],'Results_of_alltheta.mat');
        Path(kk,1)=name;
        Name(kk,1)=fullfile(seg{ii},['fish' num2str(jj)]);
        kk=kk+1;
    end
end
load(Path(15));name=fieldnames(Results_of_alltheta);
Data=struct;CR_trial={};V_trial_CS={};CR_loc={};V_trial_base={};V_trial_increased={};V_trial_preCS={};
for ii=1:length(name)
    Data=setfield(Data,name{ii},[]);
end
for zz=1:length(Path)
    if exist(Path(zz))
        disp(Path(zz))
        env=[];Results_of_alltheta=[];alltheta=[];
        load(Path(zz));ind=find(strcmp(Path,Path(zz)));
        if length(alltheta)<trial.total*frameb.per_cycle
            warning(['frist 6*fs time point was cut' ]);continue;
            %alltheta(1:(length(alltheta)/frameb.per_cycle-trial.total)*frameb.per_cycle)=[];
        end
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
            alltheta_n=normalize(alltheta,'zscore');%alltheta_n=normalize(alltheta,'range');
            dif=diff([alltheta_n(1);alltheta_n]);
            %% Vbase_spon
            env_loc_spon=[];V_base=[];
            for ii=1:length(T_spon)
                [~,env_loc,env_end,~]=env_detect(alltheta(T_spon{ii}),3*60,1);
                %h=figure;plot(alltheta(T_spon{ii}));hold on;scatter(env.env_loc(:,1),alltheta(T_spon{ii}(env.env_loc(:,1))),'filled');
                ind=T_spon{ii}([env_loc,env_end]);
                V_base_a=[];
                if ~isempty(ind)
                    for jj=1:size(ind,1)
                        V_base_a= [V_base_a;sum(abs(dif(ind(jj,1):ind(jj,2))))];
                    end
                else
                    V_base_a=[0 0];
                end
                V_base(:,ii)=sum(V_base_a)/(length(T_spon{ii})/fs.behavior);
            end
            V_base=kron(V_base(1:trial.acq_block_num+1)',ones(trial.acq_block_trial,1));
            V_base=[repmat(V_base(1,:),trial.hab(1),1);V_base;repmat(V_base(end,:),trial.hab(1)-trial.acq_block_trial,1)];
            V_trial_base{hh}(:,zz)=V_base;
            %% Vbase_preCS
            [C]=mod(env.env_loc(:,1),frameb.per_cycle);
            ind=find(C<frameb.cs_start);
            env_loc_pre_cs=[env.env_loc(ind,:),env.env_end(ind,:)];
%             figure,plot(alltheta(T_non_spon));hold on;plot(stimCS(T_non_spon)*180);
%             %scatter(env_loc_cs(:,1),alltheta(T_non_spon(env_loc_cs(:,1))),'r','filled');hold on
%             scatter(env_loc_pre_cs(:,1),alltheta(T_non_spon(env_loc_pre_cs(:,1))),'k','filled');hold on
            V_base=[];
            for ii=1:length(ind)
                V_base= [V_base;sum(abs(dif(T_non_spon(env_loc_pre_cs(ii,1)):T_non_spon(env_loc_pre_cs(ii,2)))))];
                %V_base= [V_base;sum(abs(dif(T_non_spon(env_loc_pre_cs(ii,1)):T_non_spon(env_loc_pre_cs(ii,2)))))./((env_loc_pre_cs(ii,2)-env_loc_pre_cs(ii,1))/fs.behavior)];
            end        
            [C,~,~]=unique(fix(env_loc_pre_cs(:,1)/frameb.per_cycle)+1);
            for ii=C'
                ind=find(fix(env_loc_pre_cs(:,1)/frameb.per_cycle)+1==ii);
                V_trial_preCS{hh}(ii,zz)=sum(V_base(ind))/((frameb.cs_start-1)/fs.behavior);
                 %V_trial_preCS{hh}(ii,zz)=mean(V_base(ind));
            end
            %% Vcs
            aa=(repmat(T_non_spon(env.env_loc(:,1)),length(ind_cs(:,1)),1)-repmat(ind_cs(:,1),1,length(env.env_loc(:,1))));
            ind=find(aa>=-0.1*60 & abs(aa)<=(frameb.cs_end-frameb.cs_start)+1-0.1*60);[~,ind]=ind2sub([size(aa,1),size(aa,2)],ind);
            env_loc_cs=[env.env_loc(ind,:),env.env_end(ind,:)];env_num=length(ind);
            V=[];
            for ii=1:length(ind)
                V= [V;sum(abs(dif(T_non_spon(env_loc_cs(ii,1)):T_non_spon(env_loc_cs(ii,2)))))];
                 %V= [V;sum(abs(dif(T_non_spon(env_loc_cs(ii,1)):T_non_spon(env_loc_cs(ii,2)))))./((env_loc_cs(ii,2)-env_loc_cs(ii,1))/fs.behavior)];
            end
            [C,~,~]=unique(fix(env_loc_cs(:,1)/frameb.per_cycle)+1);
            for ii=C'
                ind=find(fix(env_loc_cs(:,1)/frameb.per_cycle)+1==ii);
                V_trial_CS{hh}(ii,zz)=sum(V(ind))/((frameb.cs_end-frameb.cs_start)/fs.behavior);
                %V_trial_CS{hh}(ii,zz)=mean(V(ind));
                CR_trial{hh}(ii,zz)=length(ind);
                CR_loc{hh}(ii,zz)=mod(T_non_spon(env_loc_cs(ind(1),1)),frameb.per_cycle);
            end
        end
    else
         disp(Path(zz));warning('no mat');
    end
end
V_trial_increased_hab=mean(V_trial_increased(1:trial.hab(1),:),1);

%% test learner
%pre_cut
Vcs=V_trial_CS{3};
Vprecs=V_trial_preCS{3};
V_trial_increased=(Vcs-Vprecs)./(Vcs+Vprecs);
Vcs(find(V_trial_increased<0))=0;

P_fisher=nan(1,length(Name)); P_ranksum_Vincresed=nan(1,length(Name)); P_ranksum_Vcscs=nan(1,length(Name));P_ranksum_Vprecscs=nan(1,length(Name));
Performace_acq_slidewin=[];Performace_acq_blocks=[];
type=all;%type=[1:length(Name)];
kk=1;
for ii=type
    disp(Name(ii))
    a=Vcs(1:trial.hab(1),ii);a(find(a>0))=1;
    b=Vcs(end-trial.test(1)+1:end,ii);b(find(b>0))=1;
    pre1=length(find(a==1));pre0=length(find(a==0));
    post1=length(find(b==1));post0=length(find(b==0));
    %fisher's test
    x = table([pre1;post1],[pre0;post0],'VariableNames',{'Go','NoGo'},'RowNames',{'Pre','Post'});
    [h,p,stats] = fishertest(x,'tail','left');
    P_fisher(kk)=p;
    %ranksum  V_trial_increased
    a=V_trial_CS{3}(:,ii);
    b=V_trial_preCS{3}(:,ii);
    V_trial_increased=(a-b)./(a+b);V_trial_increased(isnan(V_trial_increased))=0;
    if sum(isnan(V_trial_increased))<length(V_trial_increased)
    [p,h,stats] = signrank(V_trial_increased(1:trial.hab(1)),V_trial_increased(end-trial.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_Vincresed(kk)=p;
    end
    %rumksum Vcs-pre VS Vcs-post
    [p,h,stats] = signrank(a(1:trial.hab(1)),a(end-trial.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_Vcscs(kk)=p;
    %rumksum Vprecs-post VS  Vcs-post 
    [p,h,stats] = signrank(b(end-trial.test(1)+1:end),a(end-trial.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_Vprecscs(kk)=p;
    %acq_slide window
    win=5;
    a=Vcs(trial.hab(1)+1:trial.hab(1)+trial.acq(1),ii);a(find(a>0))=1;b=[];
    for jj=1:trial.acq(1)-win
        b(jj)=sum(a(jj:jj+win))/win;
    end
    max(b)
    Performace_acq_slidewin(:,kk)=b;
    %acq_block
    b=[];
    for jj=1:trial.acq_block_num
        b(jj)=sum(a(1+(jj-1)*trial.acq_block_trial:(jj)*trial.acq_block_trial))/trial.acq_block_trial;
    end
    max(b)
    Performace_acq_blocks(:,kk)=b;
    kk=kk+1;
end
Performace_acq_slidewin_max=max(Performace_acq_slidewin);
Performace_acq_blocks_max=max(Performace_acq_blocks);

figure,plot(Performace_acq_slidewin(:,3))
figure,plot(Performace_acq_blocks)

%% plot
learned=[10 26 25 ];
notsure=[];
acqlearned=[ 6 7 17 8 5 18 19 23];
nonlearned=[20 3 12 13 11 14 24];
cr=[4 9 15 16 21 27];
all=[learned,acqlearned,nonlearned,cr];%[learned,notsure,acqlearned,nonlearned];
type=all;

ii=type;n=Name(ii);
a=V_trial_CS{3};a=a(:,ii);n=Name(ii);%a(end+1:42,:)=0;
b=V_trial_preCS{3};b=b(:,ii);
%V_trial_increased=a./b;
V_trial_increased=(a-b)./(a+b);
V_trial_increased(find(V_trial_increased<1/3))=nan;
figure,histogram(V_trial_increased,[-1:0.1:1]);set(gca,'fontsize',14);

Vcs=V_trial_CS{3};
Vprecs=V_trial_preCS{3};
V_trial_increased=(Vcs-Vprecs)./(Vcs+Vprecs);
b=CR_loc{3};b(find(V_trial_increased<0))=0;b=b(:,type);
b(find(b==0))=nan;b=(b-frameb.cs_start)/fs.behavior;

%state
unchanged=[6,7,8,10,15,17,18,19,20,22];
P=[23];
TP=[3,12,21,24,28];
type=[unchanged,P,TP];
b=Data.M_num_spon(1:9,type);
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

aa=a;aa(find(a<=0))=0;aa(find(a>0))=1;