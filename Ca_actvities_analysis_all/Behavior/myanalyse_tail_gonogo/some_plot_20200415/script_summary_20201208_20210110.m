%
clc;clear all;

global fs
global frameb
global trial

[fs,time,~,frameb,trial_r,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial_r.acq_block_interval=10;
trial_r.test(2)=trial_r.hab(3)+trial_r.acq_block_num*trial_r.acq_block_trial+(trial_r.acq_block_num)*trial_r.acq_block_interval+1;
trial_r.test(3)=trial_r.test(2)+trial_r.test(1)-1;
trial_r.spon_aft=[trial_r.spon_aft(1) ,trial_r.test(3)+1,trial_r.test(3)+trial_r.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial_r.acq_block_num)*fs.behavior;
trial_r.total=trial_r.total+trial_r.acq_block_interval*(trial_r.acq_block_num)
T_spon={};T_US={};l_spon={};l_US={};T_non_spon=[];
for ii=1:trial_r.acq_block_num+1
    switch ii
        case 1
            t=(trial_r.spon_bef(2)-1)*frameb.per_cycle+1:(trial_r.spon_bef(3))*frameb.per_cycle;
            l_spon{ii}='Bef.';
        case mat2cell([2:trial_r.acq_block_num],1,ones(1,trial_r.acq_block_num-1))
            t=(trial_r.hab(3)+(ii-1)*trial_r.acq_block_trial+(ii-2)*(trial_r.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial_r.hab(3)+(ii-1)*trial_r.acq_block_trial+(ii-1)*(trial_r.acq_block_interval))*frameb.per_cycle ;
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)];
        case trial_r.acq_block_num+1
            t=(trial_r.hab(3)+(ii-1)*trial_r.acq_block_trial+(ii-2)*(trial_r.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial_r.hab(3)+(ii-1)*trial_r.acq_block_trial+(ii-1)*(trial_r.acq_block_interval))*frameb.per_cycle ;
            %l_spon{ii}='Test';
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)];
        case trial_r.acq_block_num+2
            t=(trial_r.test(3))*frameb.per_cycle+1:(trial_r.total)*frameb.per_cycle;
            l_spon{ii}='Aft.';
    end
    T_spon{ii}=t;
end
% us
for ii=1:trial_r.acq_block_num+2
    switch ii
        case 1
            t=(trial_r.hab(2)-1)*frameb.per_cycle+1:(trial_r.hab(3))*frameb.per_cycle;
            l_US{ii}='Hab';
        case mat2cell([2:trial_r.acq_block_num+1],1,ones(1,trial_r.acq_block_num))
            t=(trial_r.hab(3)+(ii-2)*trial_r.acq_block_trial+(ii-2)*(trial_r.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial_r.hab(3)+(ii-1)*trial_r.acq_block_trial+(ii-2)*(trial_r.acq_block_interval))*frameb.per_cycle ;
            
            l_US{ii}=['Acq.' num2str(ii-1)];
        case trial_r.acq_block_num+2
            t=(trial_r.test(2)-1)*frameb.per_cycle+1:(trial_r.test(3))*frameb.per_cycle;
            l_US{ii}='Tst';
    end
    T_US{ii}=t;
    T_non_spon=[T_non_spon,t];
end
stimUS=zeros(1,frameb.per_cycle*trial_r.total);ind_US=[];kk=1;
for ss=1:trial_r.acq_block_num
    for tt=1:trial_r.acq_block_trial
        ind=(trial_r.hab(3)+(tt-1)+(ss-1)*(trial_r.acq_block_interval+trial_r.acq_block_trial))*frameb.per_cycle+frameb.us_start;
        stimUS(ind:ind+max(round(frameb.us_dur(ss)),1)-1)=1;%23
        ind_US(kk,1)=ind;kk=kk+1;
    end
end
trial_r_ind_CS=[trial_r.hab(2):trial_r.hab(3)];
for ii=1:trial_r.acq_block_num
    trial_r_ind_CS=[trial_r_ind_CS,[trial_r.hab(3)+1:trial_r.hab(3)+trial_r.acq_block_trial]+(ii-1)*(trial_r.acq_block_interval+trial_r.acq_block_trial)];
end
trial_r_ind_CS=[trial_r_ind_CS,trial_r.test(2):trial_r.test(3)];
stimCS=zeros(1,frameb.per_cycle*trial_r.total);ind_CS=[];kk=1;
for tt=trial_r_ind_CS
    stimCS((tt-1)*frameb.per_cycle+frameb.cs_start:(tt-1)*frameb.per_cycle+frameb.cs_end)=1;%23
    ind_CS(kk,1)=(tt-1)*frameb.per_cycle+frameb.cs_start;
    ind_CS(kk,2)=(tt-1)*frameb.per_cycle+frameb.cs_end;kk=kk+1;
end
ind_CS_hab=ind_CS(1:6,:);ind_CS_tst=ind_CS(end-5:end,:);

%load
path=uigetdir('C:\Users\zhulin\Desktop\ION\du lab\Analysing\Tail movement check\fear contioning data\');
% seg={'20201223','20201224','20201229','20201230','20201231','20210101','20210102','20210104','20210106','20210107','20210108','20210109','20210110',...
%     '20210113','20210114','20210118','20210121','20210122','20210123','20210125','20210126','20210127'};
% n=[3,4,3,3,1,1,2,1,2,2,2,2,2,1,2,2,3,1,1,1,2,1];kk=1; %(20201223-20210127,21dpf)

%seg={'20201223','20201224','20201229','20201230','20201231','20210101','20210102','20210104','20210106','20210107','20210108','20210109','20210110'};
%n=[3,4,3,3,1,1,2,1,2,2,2,2,2];kk=1; (20201223-20210110)

% seg={'20200624','20200625','20200626','20200627','20200707','20200708','20200709','20200710','20200714','20200715','20200716','20200717',...
%     '20200831','20200901','20200902','20200903','20200908','20200909','20200910','20200923','20200924','20200925','20200928','20200930',...
%     '20201003','20201019','20201020','20201021','20201025','20201026','20201027',...
%     '20201109','20201110','20201111','20201126','20201127','20201128','20201208','20201209','20201214','20201221'...
%     '20201223','20201224','20201229','20201230','20201231','20210101','20210102','20210104','20210106','20210107','20210108','20210109','20210110',...
%     '20210113','20210114','20210118','20210121','20210122','20210123','20210125','20210126','20210127'};
% n=[5,4,6,5,3,5,5,4,5,5,4,5,3,6,3,5,5,6,6,6,4,4,3,4,3,3,3,4,4,4,3,4,5,6,4,5,3,3,2,3,1,3,4,3,3,1,1,2,1,2,2,2,2,2,1,2,2,3,1,1,1,2,1];kk=1;%(20200624-20210127)

% seg={'20210223','20210224','20210309','20210311','20210401','20210402'};
% n=[3,1,1,1,2,1];kk=1;

seg={'20200626','20200627','20200909','20200928','20200930','20201027','20201229','20210107','20210108','20210127','20210402'};
n=[4,5,4,2,1,3,1,1,2,1,1];kk=1;


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
Data=struct;CR_trial_r={};CR_precs_trial_r={};
V_trial_r_CS={};CR_loc={};V_trial_r_base={};V_trial_r_increased={};V_trial_r_preCS={};CR_trial_r_ind={};
env_loc_non_spon={};
P_bino_movratio_cscs_b=[];P_bino_movratio_precscs_post_b=[];P_bino_movratio_sponspon_b=[];P_bino_movratio_cscs_dur_b=[];P_bino_movratio_precscs_dur_b=[];
P_bino_movratio_cscs=[];P_bino_movratio_precscs_post=[];P_bino_movratio_sponspon=[]; P_bino_movratio_cscs_dur=[];P_bino_movratio_precscs_dur=[];P_bino_movratio_cscs_dur_min=[];P_bino_movratio_precscs_dur_min=[];
for ii=1:length(name)
    Data=setfield(Data,name{ii},[]);
end
for zz=1:length(Path)
    if exist(Path(zz))
        disp(Path(zz))
        env=[];Results_of_alltheta=[];alltheta=[];
        load(Path(zz));ind=find(strcmp(Path,Path(zz)));
        if length(alltheta)<trial_r.total*frameb.per_cycle
            warning(['less frame' ]);continue;
            %alltheta(1:(length(alltheta)/frameb.per_cycle-trial_r.total)*frameb.per_cycle)=[];
        end
         if length(alltheta)>trial_r.total*frameb.per_cycle
            warning(['more frame' ]);continue;
            %alltheta(1:(length(alltheta)/frameb.per_cycle-trial_r.total)*frameb.per_cycle)=[];
        end
        for ii=1:length(name)
            if min(size(getfield(Results_of_alltheta,name{ii})))<=1
                b=getfield(Data,name{ii})';a=getfield(Results_of_alltheta,name{ii});if ~isempty(a) b(ind,1:length(a))=a;end
            elseif min(size(getfield(Results_of_alltheta,name{ii})))>1
                b=getfield(Data,name{ii})';
                for jj=1:trial_r.acq_block_trial
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
        %%
        startpoint=env.env_loc;re_startpoint=[];env_loc_non_spon{zz}=startpoint;
        re_startpoint(:,1)=fix(startpoint/frameb.per_cycle)+1;
        re_startpoint(:,2)=mod(startpoint,frameb.per_cycle);
        a=0;trial.spon_bef=[a min(1,a) a];
        a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
        trial.acq_block_num=8;
        trial.acq_block_trial=3;
        trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
        a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
        a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
        trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);
        isbinary=1;
        [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(alltheta(T_non_spon)',re_startpoint,[],[],[],isbinary);
        % pre_con=[bef_cond.befCS_shake/2 bef_cond.CS_shake];
        % aft_con=[aft_cond.befCS_shake/2 aft_cond.CS_shake];
        P_bino_movratio_cscs_b(:,zz)=p_value.CS_15_6;
        P_bino_movratio_precscs_post_b(:,zz)=p_value.aftcond_befCS_CS_6;
        P_bino_movratio_sponspon_b(:,zz)=p_value.states_6;
        [P_bino_movratio_cscs_dur_b_min(:,zz),I]=min(p_value.dur_learnning_CS_6);
        P_bino_movratio_precscs_dur_b_min(:,zz)=p_value.durcond_befCS_CS(I);
        P_bino_movratio_cscs_dur_b(:,zz)=p_value.dur_learnning_CS_6;
        P_bino_movratio_precscs_dur_b(:,zz)=p_value.durcond_befCS_CS;
        isbinary=0;
        [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(alltheta(T_non_spon)',re_startpoint,[],[],[],isbinary);
        P_bino_movratio_cscs(:,zz)=p_value.CS_15_6;
        P_bino_movratio_precscs_post(:,zz)=p_value.aftcond_befCS_CS_6;
        P_bino_movratio_sponspon(:,zz)=p_value.states_6;
        [P_bino_movratio_cscs_dur_min(:,zz),I]=min(p_value.dur_learnning_CS_6);
        P_bino_movratio_precscs_dur_min(:,zz)=p_value.durcond_befCS_CS(I);
        P_bino_movratio_cscs_dur(:,zz)=p_value.dur_learnning_CS_6;
        P_bino_movratio_precscs_dur(:,zz)=p_value.durcond_befCS_CS;

        %%
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
            V_base=kron(V_base(1:trial_r.acq_block_num+1)',ones(trial_r.acq_block_trial,1));
            V_base=[repmat(V_base(1,:),trial_r.hab(1),1);V_base;repmat(V_base(end,:),trial_r.hab(1)-trial_r.acq_block_trial,1)];
            V_trial_r_base{hh}(:,zz)=V_base;
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
                V_trial_r_preCS{hh}(ii,zz)=sum(V_base(ind))/((frameb.cs_start-1)/fs.behavior);
                %V_trial_r_preCS{hh}(ii,zz)=mean(V_base(ind));
                CR_precs_trial_r{hh}(ii,zz)=length(ind);
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
                V_trial_r_CS{hh}(ii,zz)=sum(V(ind))/((frameb.cs_end-frameb.cs_start)/fs.behavior);
                %V_trial_r_CS{hh}(ii,zz)=mean(V(ind));
                CR_trial_r{hh}(ii,zz)=length(ind);
                CR_trial_r_ind{hh}{ii,zz}=mod(T_non_spon(env_loc_cs(ind,1)),frameb.per_cycle);
                CR_loc{hh}(ii,zz)=mod(T_non_spon(env_loc_cs(ind(1),1)),frameb.per_cycle);
            end
        end
    else
        disp(Path(zz));warning('no mat');
    end
end

%% P value of wilcoxon & fisher of V
P_fisher=nan(1,length(Name)); P_ranksum_Vincresed=nan(1,length(Name)); P_ranksum_Vcscs=nan(1,length(Name));P_ranksum_Vprecscs_post=nan(1,length(Name));P_ranksum_Vsponspon=nan(1,length(Name));
Performace_acq_slidewin=[];Performace_acq_blocks=[];P_ranksum_Vprecscs_dur=[];P_ranksum_Vcscs_dur=[];
Performace_hab_slidewin=[];Performace_test_slidewin=[];
kk=1;
for ii=1:size(Name,1)
    disp(Name(ii))
    
    %fisher's test
    Vcs=V_trial_r_CS{3};
    a=Vcs(1:trial_r.hab(1),ii);a(find(a>0))=1;
    b=Vcs(end-trial_r.test(1)+1:end,ii);b(find(b>0))=1;
    pre1=length(find(a==1));pre0=length(find(a==0));
    post1=length(find(b==1));post0=length(find(b==0));
    x = table([pre1;post1],[pre0;post0],'VariableNames',{'Go','NoGo'},'RowNames',{'Pre','Post'});
    [h,p,stats] = fishertest(x,'tail','left');
    P_fisher(kk)=p;
    
    %ranksum  V_trial_r_increased
    a=V_trial_r_CS{3}(:,ii);
    b=V_trial_r_preCS{3}(:,ii);
    V_trial_r_increased=(a-b)./(a+b);V_trial_r_increased(isnan(V_trial_r_increased))=0;
    if sum(isnan(V_trial_r_increased))<length(V_trial_r_increased)
        [p,h,stats] = signrank(V_trial_r_increased(1:trial_r.hab(1)),V_trial_r_increased(end-trial_r.test(1)+1:end),'alpha',0.05,'tail','left');
        P_ranksum_Vincresed(kk)=p;
    end
    
    %rumksum Vcs-pre VS Vcs-post
    [p,h,stats] = ranksum(a(1:trial_r.hab(1)),a(end-trial_r.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_Vcscs(kk)=p;
    %rumksum Vprecs-post VS  Vcs-post
    [p,h,stats] = signrank(b(end-trial_r.test(1)+1:end),a(end-trial_r.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_Vprecscs_post(kk)=p;
    %rumksum Vprecs-pre VS Vprecs-post
    [p,h,stats] = ranksum(b(1:trial_r.hab(1)),b(end-trial_r.test(1)+1:end),'alpha',0.05);
    P_ranksum_Vsponspon(kk)=p;
    % rumksum Vcs-pre VS Vcs-dur  
    win=6;
    for zz=trial.acq(2):trial.acq(3)-win
        [p,h,stats] = ranksum(a(1:trial.hab(1)),a(zz:zz+win-1),'alpha',0.05,'tail','left');
        P_ranksum_Vcscs_dur(zz,kk)=p;
    end
    %rumksum Vprecs-dur VS  Vcs-dur
    for zz=trial.acq(2):trial.acq(3)-win
        [p,h,stats] = signrank(b(zz:zz+win-1),a(zz:zz+win-1),'alpha',0.05,'tail','left');
        P_ranksum_Vprecscs_dur(zz,kk)=p;
    end
    
    %Performace acq_slide window
    win=6;
    a=Vcs(trial_r.hab(1)+1:trial_r.hab(1)+trial_r.acq(1),ii);a(find(a>0))=1;b=[];
    for jj=1:trial_r.acq(1)-win
        b(jj)=sum(a(jj:jj+win-1))/win;
    end
    max(b)
    Performace_acq_slidewin(:,kk)=b;
    %acq_block
    b=[];
    for jj=1:trial_r.acq_block_num
        b(jj)=sum(a(1+(jj-1)*trial_r.acq_block_trial:(jj)*trial_r.acq_block_trial))/trial_r.acq_block_trial;
    end
    max(b)
    Performace_acq_blocks(:,kk)=b;
    
    a=Vcs(trial.hab(2):trial.hab(3),ii);a(find(a>0))=1;
    Performace_hab_slidewin(:,kk)=sum(a)/trial_r.hab(1);
    
    a=Vcs(trial.test(2):trial.test(3),ii);a(find(a>0))=1;
    Performace_test_slidewin(:,kk)=sum(a)/trial_r.hab(1);
    
    kk=kk+1;
end
Performace_acq_slidewin_max=max(Performace_acq_slidewin);
Performace_acq_blocks_max=max(Performace_acq_blocks);

figure,plot(Performace_acq_slidewin(:,3))
figure,plot(Performace_acq_blocks)

%% P value of wilcoxon & fisher of CrNUM
P_ranksum_CRnum=nan(1,length(Name)); P_ranksum_CRnumcscs=nan(1,length(Name));P_ranksum_CRnumprecscs_post=nan(1,length(Name));P_ranksum_CRnumsponspon=nan(1,length(Name));
P_ranksum_CRnumprecscs_dur=[];P_ranksum_CRnumcscs_dur=[];
kk=1;
for ii=1:size(Name,1)
    disp(Name(ii))
    b=CR_precs_trial_r{3}(:,ii);
    a=CR_trial_r{3}(:,ii);
    %rumksum CRnum cs-pre VS CRnum cs-post
    [p,h,stats] = ranksum(a(1:trial_r.hab(1)),a(end-trial_r.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_CRnumcscs(kk)=p;
    %rumksum CRnum precs-post VS  CRnum cs-post
    [p,h,stats] = signrank(b(end-trial_r.test(1)+1:end),a(end-trial_r.test(1)+1:end),'alpha',0.05,'tail','left');
    P_ranksum_CRnumprecscs_post(kk)=p;
    %rumksum CRnum precs-pre VS CRnum precs-post
    [p,h,stats] = ranksum(b(1:trial_r.hab(1)),b(end-trial_r.test(1)+1:end),'alpha',0.05);
    P_ranksum_CRnumsponspon(kk)=p;
    % rumksum CRnum cs-pre VS CRnum cs-dur  
    win=6;
    for zz=trial.acq(2):trial.acq(3)-win
        [p,h,stats] = ranksum(a(1:trial.hab(1)),a(zz:zz+win-1),'alpha',0.05,'tail','left');
        P_ranksum_CRnumcscs_dur(zz,kk)=p;
    end
    %rumksum CRnum precs-dur VS  CRnum cs-dur
    for zz=trial.acq(2):trial.acq(3)-win
        [p,h,stats] = signrank(b(zz:zz+win-1),a(zz:zz+win-1),'alpha',0.05,'tail','left');
        P_ranksum_CRnumprecscs_dur(zz,kk)=p;
    end
    
    kk=kk+1;
end

%% plot

%20200624-20210127
learned=[13,17,20,107,109,83,134,135,179,186,192,195];%12
acqlearned=[57,89,33,78,161,160,153,154,145,25,174,194,202,204];%14
nonlearned=[75,77,87,88,106,112,73,7,9,28,5,18,122,123,121,141,142,143,140,158,138,147,148,146,149,159,139,150,...
             177,181,189,190,191,193,199,203,206,208,209,211,212 ];%41
confusing=[213];

%20201223-20210127 
% learned=[8,15,21,24];
% acqlearned=[3,31,33];
% nonlearned=[6,7,10,12,17,18,19,20,22,23,28,32,35,37,38,40,41,42]; 
% cr=[];

% all=[learned,acqlearned,nonlearned,confusing];[1:length(Name)];%
all=[4,6,9,13,15,16,18,19,20,21,23,24,25]
type=all;%type=[1:length(Name)];
n=Name(type);
P_fisher_index=P_fisher(:,type);

Performace_acq_slidewin_max_index=Performace_acq_slidewin_max(:,type);
Performace_acq_blocks_max_index=Performace_acq_blocks_max(:,type);
Performace_acq_slidewin_index=Performace_acq_slidewin(:,type);%not max. itis learning performance

P_bino_movratio_cscs_index=P_bino_movratio_cscs(:,type);
P_bino_movratio_precscs_post_index=P_bino_movratio_precscs_post(:,type);
P_bino_movratio_sponspon_index=P_bino_movratio_sponspon(:,type);
P_bino_movratio_cscs_dur_index=P_bino_movratio_cscs_dur(:,type);
P_bino_movratio_precscs_dur_index= P_bino_movratio_precscs_dur(:,type);

P_bino_movratio_cscs_b_index=P_bino_movratio_cscs_b(:,type);
P_bino_movratio_precscs_post_b_index=P_bino_movratio_precscs_post_b(:,type);
P_bino_movratio_sponspon_b_index=P_bino_movratio_sponspon_b(:,type);
P_bino_movratio_cscs_dur_b_index=P_bino_movratio_cscs_dur_b(:,type);
P_bino_movratio_precscs_dur_b_index= P_bino_movratio_precscs_dur_b(:,type);

P_ranksum_Vincresed_index=P_ranksum_Vincresed(:,type);
P_ranksum_Vcscs_index=P_ranksum_Vcscs(:,type);
P_ranksum_Vprecscs_post_index=P_ranksum_Vprecscs_post(:,type);
P_ranksum_Vsponspon_index=P_ranksum_Vsponspon(:,type);
%P_ranksum_Vcscs_dur and P_ranksum_Vprecscs_dur in line 293

ii=type;n=Name(ii);
a=V_trial_r_CS{3};a=a(:,ii);n=Name(ii);%a(end+1:42,:)=0;
b=V_trial_r_preCS{3};b=b(:,ii);
%V_trial_r_increased=a./b;
V_trial_r_increased=(a-b)./(a+b);
%V_trial_r_increased(find(V_trial_r_increased<1/3))=nan;
figure,histogram(V_trial_r_increased,[-1:0.1:1]);set(gca,'fontsize',14);

% Vcs=V_trial_r_CS{3};
% Vprecs=V_trial_r_preCS{3};
% V_trial_r_increased=(Vcs-Vprecs)./(Vcs+Vprecs);
b=CR_loc{3};
%b(find(V_trial_r_increased<0))=0;
b=b(:,learned);
b(find(b==0))=inf;b=(b-frameb.cs_start)/fs.behavior;b=flip(b);

%state
%20200623-20210127
unchanged=[177,178,181,186,189,19,20,22,33,37,38];%11
Active=[];
P=[23,40,41];
TP=[3,21,24,28];
TA=[31,32,35];

% %20201223-20210127
% unchanged=[6,8,10,15,18,19,20,22,33,37,38];%11
% Active=[];
% P=[23,40,41];
% TP=[3,21,24,28];
% TA=[31,32,35];

type=[unchanged,Active,P,TP,TA];
b=Data.M_num_spon(1:9,type);
bb=(b-repmat(b(1,:),9,1))./(b+repmat(b(1,:),9,1));

a=V_trial_r_base{3};a=reshape(a,1,[]);
b=V_trial_r_preCS{3};b=reshape(b,1,[]);
c=V_trial_r_CS{3};c=reshape(c,1,[]);
figure,plot([a;b;c],'k--.');xlim([0.5 3.5])

%% mov. prob
bin_win=0.4*fs.behavior; bin=reshape(1:frameb.per_cycle,bin_win,[]);
num={};
for zz=1:3
    re_startpoint_1=[];re_startpoint_2=[];
    switch zz
        case 1
            type=learned;
        case 2
            type=acqlearned;
        case 3
            type=nonlearned;
    end
    for ii=type
        startpoint=env_loc_non_spon{ii};
        re_startpoint_1=[re_startpoint_1;fix(startpoint/frameb.per_cycle)+1];
        re_startpoint_2=[re_startpoint_2;mod(startpoint,frameb.per_cycle)];
    end
    re_startpoint=[re_startpoint_1,re_startpoint_2];
    for ii=1:trial.total
        ind=find(re_startpoint(:,1)==ii);
        for jj=1:size(bin,2)
            num{zz}(ii,jj)=length(intersect(re_startpoint(ind,2),bin(:,jj)))/length(type);
        end
    end
end
t=([1:bin_win:frameb.per_cycle]-1)/fs.behavior;
t(find(t~=12 & t~= 16.8 & t~=30 & t~=0))=nan;
for zz=1:3
    a=num{zz};
    a=[a(1:trial.hab(1),:);nan(1,size(a,2));a(trial.acq(2):trial.acq(3),:);nan(1,size(a,2));a(trial.test(2):trial.test(3),:)];
    num{zz}=a;
end