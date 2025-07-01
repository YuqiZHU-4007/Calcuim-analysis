clc;
clear all;
labels_typei_ADJ={'CSavtivation','CSinhibition','USactivation','CSupregulate','CSdownregulate','CSstableregulate','USupregulate','USdownregulate','USstableregulate','CSnewemergedneuron'};
labels_typei_other={'NUM_UP_UR','NUM_UP_CR','NUM_DOWN_CR','CS_up','CS_down','CS_stable','US_up','US_down','US_stable','NUM_NEW_EMERGED'};
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'});
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath,'/Path.mat'])
iscutmov=1;
res=[0.66,0.66,10];
load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);

A=struct;
A=setfield(A,'NUM_UP_UR',[]);A=setfield(A,'NUM_UP_CR',[]);A=setfield(A,'NUM_DOWN_CR',[]);
A=setfield(A,'CS_newemergedneuron',[]);A=setfield(A,'NUM_NEW_EMERGED',[]);
A=setfield(A,'CS_up',[]);A=setfield(A,'CS_down',[]);A=setfield(A,'CS_stable',[]);
A=setfield(A,'US_up',[]);A=setfield(A,'US_down',[]);A=setfield(A,'US_stable',[]);
A=setfield(A,'CS_newemergedneuron',[]);
A=setfield(A,'ind_all_CSUS_RESPONSIVE',[]);
% for type=1:length(labels_typei_ADJ)
% A=setfield(A,labels_typei_ADJ{type},{});
% end
path_all={};kk=1;
for ii=1:4
    path=Path{ii};
    for jj=1:length(path)
        p=path{jj};path_all{kk}=p;kk=kk+1;
        load(fullfile(p,'CR_ind_summary.mat'));
        a=getfield(A,'NUM_UP_UR');a=cat(2,a,NUM_UP_UR'./size(CR_ind_up,1));A=setfield(A,'NUM_UP_UR',a);
        a=getfield(A,'NUM_UP_CR');a=cat(2,a,NUM_UP_CR(:,iscutmov)./size(CR_ind_up,1));A=setfield(A,'NUM_UP_CR',a);
        a=getfield(A,'NUM_DOWN_CR');a=cat(2,a,NUM_DOWN_CR(:,iscutmov)./size(CR_ind_up,1));A=setfield(A,'NUM_DOWN_CR',a);
        a=getfield(A,'CS_up');a=cat(2,a,CS_up(:,iscutmov));A=setfield(A,'CS_up',a);
        a=getfield(A,'CS_down');a=cat(2,a,CS_down(:,iscutmov));A=setfield(A,'CS_down',a);
        a=getfield(A,'CS_stable');a=cat(2,a,CS_stable(:,iscutmov));A=setfield(A,'CS_stable',a);
        a=getfield(A,'US_up');a=cat(2,a,US_up(:,iscutmov));A=setfield(A,'US_up',a);
        a=getfield(A,'US_down');a=cat(2,a,US_down(:,iscutmov));A=setfield(A,'US_down',a);
        a=getfield(A,'US_stable');a=cat(2,a,US_stable(:,iscutmov));A=setfield(A,'US_stable',a);
        for fraction_type=1:2
            for type=1:length(labels_typei_ADJ)
                f=['Fraction_of_ind_all_CSUS_RESPONSIVE',num2str(fraction_type),labels_typei_ADJ{type}];
                if isfield(A,f)
                    a=getfield(A,f);
                    a=cat(3,a,Fraction_in_region_type{type,iscutmov,fraction_type});
                    A=setfield(A,f,a);
                else
                    A=setfield(A,f,{});
                end
            end
        end
        for fraction_type=1:2
            for type=1:length(labels_typei_ADJ)
                for type_ind=1:3
                    f=['Fraction_of_ind_all_CSUS_RESPONSIVE',num2str(fraction_type),labels_typei_ADJ{type}];
                    if isfiled(A,f)
                        a=getfield(A,f);
                        a=cat(3,a,Fraction_in_region_type_emerged{type,iscutmov,type_ind,fraction_type});
                        A=setfield(A,f,a);
                    else
                        A=setfield(f,{});
                    end
                end
            end
        end
        for type=1:length(labels_typei_ADJ)
            aa={};
            for session=1:size(ind_all_CSUS_RESPONSIVE,2)
                aa=cat(1,aa,ind_all_CSUS_RESPONSIVE{type,session,iscutmov});
            end
            f=['ind_all_CSUS_RESPONSIVE',labels_typei_ADJ{type}];
            if isfield(A,f)
                a=getfield(A,f);a=cat(2,a,aa);
                A=setfield(A,f,a);
            else
                A=setfield(A,f,{});
            end
        end
        for fraction_type=1:3
            for type=1:length(labels_typei_ADJ)
                aa={};
                for session=1:size(ind_all_emergeCSUS,2)
                    aa=cat(1,aa,ind_all_emergeCSUS{type,session,iscutmov,fraction_type});
                end
                f=['ind_all_emergeCSUS',labels_typei_ADJ{type},num2str(fraction_type)];
                if isfield(A,f)
                    a=getfield(A,f);a=cat(2,a,aa);
                    A=setfield(A,f,a);
                else
                    A=setfield(A,f,[]);
                end
            end
        end
    end
end

%         nn=[p(end-14:end-7),p(end-5:end-1)];
%         supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
%         ind=find(CR_ind_up(:,1,iscutmov)==1);
%         ind2=find(CR_ind_down(:,1,iscutmov)==1);
%         ind3=union(ind,ind2);
%         ind3=setdiff(1:size(CR_ind_up(:,1,iscutmov),1),ind3);Fraction_in_region_type=[];NUM_ind=[];
%         for zz=2:length(sessionx)
%             switch zz
%                 case 10
%                     ind2=find(CR_ind_up(:,2,iscutmov)==1);
%                 otherwise
%                     ind2=find(CR_ind_up_acq(:,zz-1,iscutmov)==1);
%             end
%             CS_new_emerged=intersect(ind2,ind3);
%             ind=CS_new_emerged;
%             NUM_ind(:,zz)=length(ind)./size(CR_ind_up,1);
%             gIX=ones(1,length(ind))';supervoxel_in_fish=supervolxeli(ind,1:3);clrmap = hsv(max(gIX));
%             [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervoxel_in_fish,reg_mask,reg_name,reg_loc,temp_supervoxel,temp_env,clrmap,nn);
%             Fraction_in_region_type(:,zz)=fraction_in_region_in_clust(:,1);
%         end
%         a=getfield(A,'NUM_NEW_EMERGED');a=cat(2,a,NUM_ind');A=setfield(A,'NUM_NEW_EMERGED',a);
%         a=getfield(A,'CS_newemergedneuron');a=cat(3,a,Fraction_in_region_type);A=setfield(A,'CS_newemergedneuron',a);


numL=[1,2,3,4,6,7];
numC=8:10;numNL=11:16;numFL=17:20;
group={'L','FL','NL','C'};
regionxx=categorical(Label_region{1,1}(:,1));regionxx = reordercats(regionxx,Label_region{1,1}(:,1));
mean1= @(x)(mean(x,2));
error1= @(x)(std(x,[],2)./sqrt(length(x)));error2= @(x)(mean(x,2)-std(x,[],2));
statistic_test=@(x,y)(ranksum(x,y));
clr=[1 0 0;1 0.41 0.16;0 0 1;0.15 0.15 0.15];
%labels_typei_ADJ={'CSavtivation','CSinhibition','USactivation','CSupregulate','CSdownregulate','CSstableregulate','USupregulate','USdownregulate','USstableregulate','CS_newemergedneuron'};
%% 3
for type=2:length(labels_typei_other)
    h=figure('position',[33,129,1721,403],'name',labels_typei_other{type},'color','w');
    tiledlayout(1,4);
    data=getfield(A,labels_typei_other{type});
    for kk=1:4
        switch kk
            case 1
                ind=numL;
            case 4
                ind=numC;
            case 3
                ind=numNL;
            case 2
                ind=numFL;
        end
        dataL=squeeze(data(:,ind));
        ax1 = nexttile;
        b=bar(sessionx,mean1(dataL));hold on;
        b.CData = clr(kk,:);
        b.FaceColor = 'flat';
        b.FaceAlpha = 0.5;
        er = errorbar(sessionx,mean1(dataL),error1(dataL));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        title([labels_typei_other{type},group(kk)],'Interpreter','none');set(gca,'fontsize',16);
        xtips1 = b.XEndPoints;ytips1 = b.YEndPoints;
        for jj=1:size(dataL,1)
            if type>3
                n=2;
            else
                n=1;
            end
            p=signrank(dataL(n,:),dataL(jj,:));
            if p<0.05
                text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(4,:),'fontsize',14);hold on;
                disp([labels_typei_other{type},group{kk},'----SESSION ',num2str(jj),'p<0.05!!!!!!']);
            end
        end
    end
    saveas(h,fullfile(savepath,strcat(string(labels_typei_other{type}),num2str(1))),'jpeg');
    close(h);
end
%% 4
for type=2:length(labels_typei_other)
    data=getfield(A,labels_typei_other{type});
    dataL=squeeze(data(:,numL));
    dataC=squeeze(data(:,numC));
    dataNL=squeeze(data(:,numNL));
    dataFL=squeeze(data(:,numFL));
    datam = cat(2,mean1(dataL),mean1(dataFL),mean1(dataNL),mean1(dataC))';
    h=figure('position',[33,129,1721,403],'name',labels_typei_other{type},'color','w');
    b=bar(sessionx,datam);title(labels_typei_other{type},'Interpreter','none');
    for k = 1:size(b,2)
        b(k).CData = clr(k,:);
        b(k).FaceColor = 'flat';
        b(k).FaceAlpha = 0.5;
    end
    set(gca,'fontsize',16)
    hold on;
    legend(group);
    xtips1 = b(1).XEndPoints;xtips2 = b(2).XEndPoints; xtips3 = b(3).XEndPoints;xtips4 = b(4).XEndPoints;
    ytips1 = b(1).YEndPoints;ytips2 = b(2).YEndPoints; ytips3 = b(3).YEndPoints;ytips4 = b(4).YEndPoints;
    for jj=1:length(sessionx)
        p=ranksum(dataL(jj,:),dataC(jj,:));
        if  p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'L-C','p<0.05!!!!!!']);
            text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','top','color',clr(2,:),'fontsize',14);hold on;
        end
        p=ranksum(dataL(jj,:),dataNL(jj,:));
        if  p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'L_NL','p<0.05!!!!!!']);
            text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','middle','color',clr(3,:),'fontsize',14);hold on;
        end
        p=ranksum(dataL(jj,:),dataFL(jj,:));
        if  p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'L_FL','p<0.05!!!!!!']);
            text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(4,:),'fontsize',14);hold on;
        end
        p=ranksum(dataC(jj,:),dataNL(jj,:));
        if p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'C_NL','p<0.05!!!!!!']);
            text(xtips3(jj),ytips3(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(2,:),'fontsize',14);hold on;
        end
        p=ranksum(dataC(jj,:),dataFL(jj,:));
        if  p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'C_FL','p<0.05!!!!!!']);
            text(xtips2(jj),ytips2(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(2,:),'fontsize',14);hold on;
        end
        p=ranksum(dataNL(jj,:),dataFL(jj,:));
        if  p<=0.05;
            disp([labels_typei_other{type},'----SESSION ',num2str(jj),'NL_FL','p<0.05!!!!!!']);
            text(xtips2(jj),ytips2(jj),'*','HorizontalAlignment','center','VerticalAlignment','middle','color',clr(3,:),'fontsize',14);hold on;
        end
    end
    saveas(h,fullfile(savepath,strcat(string(labels_typei_other{type}),num2str(2))),'jpeg');
    close(h);
end

%% 1
pL_C={};pL_NL={};pL_FL={};pC_NL={};pC_FL={};pNL_FL={};data_ALL={};
for type=1:length(labels_typei_ADJ)
    dataall=getfield(A,labels_typei_ADJ{type});
    for ii=1:size(dataall,2)
        dataL=squeeze(dataall(:,ii,numL));
        dataC=squeeze(dataall(:,ii,numC)); 
        %dataC(5:6,2:3)=dataC(5:6,2:3)+0.02;
        dataNL=squeeze(dataall(:,ii,numNL));
        dataFL=squeeze(dataall(:,ii,numFL));
        ll=max([size(dataL,2),size(dataFL,2),size(dataNL,2),size(dataC,2)]);
        data_ALL{type,ii} = cat(2,dataL,[dataFL,nan(size(dataFL,1),ll-size(dataFL,2))],[dataNL,nan(size(dataNL,1),ll-size(dataNL,2))],[dataC,nan(size(dataC,1),ll-size(dataC,2))]);
        %data_ALL{type,ii}(5:6,20:21)=data_ALL{type,ii}(5:6,20:21)+0.15;
        data = cat(2,mean1(dataL),mean1(dataFL),mean1(dataNL),mean1(dataC))';
        errhigh = cat(2,error1(dataL),error1(dataFL),error1(dataNL),error1(dataC))';
        %subplot(5,2,ii);
        h=figure('position',[33,129,1721,403],'name',labels_typei_ADJ{type},'color','w');
        b=bar(regionxx,data);title([labels_typei_ADJ{type},sessionx(ii)],'Interpreter','none');
        for k = 1:size(b,2)
            b(k).CData = clr(k,:);
            b(k).FaceColor = 'flat';
            b(k).FaceAlpha = 0.5;
        end
        set(gca,'fontsize',16)
        %ylim([0 0.1])
        hold on;
        er = errorbar(data,errhigh);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        legend(group);
        xtips1 = b(1).XEndPoints;xtips2 = b(2).XEndPoints; xtips3 = b(3).XEndPoints;xtips4 = b(4).XEndPoints;
        ytips1 = b(1).YEndPoints;ytips2 = b(2).YEndPoints; ytips3 = b(3).YEndPoints;ytips4 = b(4).YEndPoints;
        for jj=1:length(regionxx)
            pL_C{type}(ii,jj)=ranksum(dataL(jj,:),dataC(jj,:));
            if  pL_C{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'L-C','p<0.05!!!!!!']);
                %ind=find(pL_C{type}(ii,jj)<0.05);
                %a=xtips1(ind);b=xtips2(ind);c=max([ytips1(ind)',ytips2(ind)']);
                %line([a,b],[c,c])
                text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','top','color',clr(2,:),'fontsize',14);hold on;
            end
            pL_NL{type}(ii,jj)=ranksum(dataL(jj,:),dataNL(jj,:));
            if  pL_NL{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'L_NL','p<0.05!!!!!!']);
                text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','middle','color',clr(3,:),'fontsize',14);hold on;
            end
            pL_FL{type}(ii,jj)=ranksum(dataL(jj,:),dataFL(jj,:));
            if  pL_FL{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'L_FL','p<0.05!!!!!!']);
                text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(4,:),'fontsize',14);hold on;
            end
            pC_NL{type}(ii,jj)=ranksum(dataC(jj,:),dataNL(jj,:));
            if pC_NL{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'C_NL','p<0.05!!!!!!']);
                text(xtips3(jj),ytips3(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(3,:),'fontsize',14);hold on;
            end
            pC_FL{type}(ii,jj)=ranksum(dataC(jj,:),dataFL(jj,:));
            if  pC_FL{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'C_FL','p<0.05!!!!!!']);
                text(xtips2(jj),ytips2(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(2,:),'fontsize',14);hold on;
            end
            pNL_FL{type}(ii,jj)=ranksum(dataNL(jj,:),dataFL(jj,:));
            if  pNL_FL{type}(ii,jj)<=0.05;
                disp([labels_typei_ADJ(type),'----SESSION ',num2str(ii),regionxx(jj),'NL_FL','p<0.05!!!!!!']);
                text(xtips2(jj),ytips2(jj),'*','HorizontalAlignment','center','VerticalAlignment','middle','color',clr(3,:),'fontsize',14);hold on;
            end
        end
        saveas(h,fullfile(savepath,strcat(string(labels_typei_ADJ{type}),'_',string(sessionx(ii)),'.jpeg')));
        close(h);
    end
end
%% 2
p_session={};
for type=1:length(labels_typei_ADJ)
    dataall=getfield(A,labels_typei_ADJ{type});
    for ii=1:size(dataall,1)
        h=figure('position',[33,129,1721,403],'name',labels_typei_ADJ{type},'color','w');
        tiledlayout(1,4);
        for kk=1:4
            switch kk
                case 1
                    ind=numL;
                case 4
                    ind=numC;
                case 3
                    ind=numNL;
                case 2
                    ind=numFL;
            end
            dataL=squeeze(dataall(ii,:,ind));%dataL=dataL+0.15;
%             if kk==4
%                 dataL=dataL+0.02;
%             end
            ax1 = nexttile;
            b=bar(sessionx(1:size(dataL,1)),mean1(dataL));hold on;
            b.CData = clr(kk,:);
            b.FaceColor = 'flat';
            b.FaceAlpha = 0.5;
            er = errorbar(sessionx(1:size(dataL,1)),mean1(dataL),error1(dataL));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            title([labels_typei_ADJ{type},regionxx(ii),group(kk)],'Interpreter','none');set(gca,'fontsize',16);
            xtips1 = b.XEndPoints;ytips1 = b.YEndPoints;
            %ylim([0 0.12])
            for jj=1:size(dataL,1)
                if type>3 && type<=6 ||  type==10
                    n=2;
                elseif type>6 && type<=9
                    n=3;
                else 
                    n=1;
                end
                p_session{type,kk}(jj,ii)=signrank(dataL(n,:),dataL(jj,:));
                if p_session{type,kk}(jj,ii)<0.05
                    text(xtips1(jj),ytips1(jj),'*','HorizontalAlignment','center','VerticalAlignment','bottom','color',clr(2,:),'fontsize',14);hold on;
                    disp([labels_typei_ADJ(type),group{kk},regionxx(ii),'----SESSION ',num2str(jj),'p<0.05!!!!!!']);
                end
            end
        end
         saveas(h,fullfile(savepath,strcat(string(labels_typei_ADJ{type}),'_',string(regionxx(ii))),''),'jpeg');
         close(h);
    end
end
%% 
explode = [1 1 0 1];label= {'Up','Down','Stable','None'};
for kk=1:4
    h=figure('position',[14,187,1833,716]);
    switch kk
        case 1
            ind=numL;
        case 4
            ind=numC;
        case 3
            ind=numNL;
        case 2
            ind=numFL;
    end
    for ii=1:10
        CS_up_l=A.CS_up(ii,ind)';CS_down_l=A.CS_down(ii,ind)';CS_stable_l=A.CS_stable(ii,ind)';CS_none_l=1-CS_up_l-CS_down_l-CS_stable_l;
        subplot(1,10,ii),
        pie(mean([CS_up_l CS_down_l CS_stable_l CS_none_l],1),explode);%title(sessionx(ii),'VerticalAlignment','bottom');
        set(gca,'fontsize',14,'FontWeight','bold','linewidth',2)
        legend(label)
    end
end
%% behavioral CR_ratio
 CR_ratio={};
for zz=1:4
    path=Path{zz};
    for jj=1:length(path)
        p=path{jj};
        load(fullfile(p,'para.mat'));
        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
        for ii=1:trial.acq_block_num+2
            switch ii
                case 1
                    trial_ind_session=trial.hab(2):trial.hab(3);
                case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
                    trial_ind_session=(trial.hab(3)+(ii-2)*trial.acq_block_trial)+1 :(trial.hab(3)+(ii-1)*trial.acq_block_trial);
                case trial.acq_block_num+2
                    trial_ind_session=trial.test(2):trial.test(3);
            end
            CR_ratio{zz}(ii,jj)=length(unique(CR(find(CR(:,1)<=trial_ind_session(end) & CR(:,1)>=trial_ind_session(1)),1)))/length( trial_ind_session);
        end
    end
end
%% ind_all_CSUS_RESPONSIVE
%ind
%location
%co-localtion
for batchi=1:4
    path=Path{batchi};loc={};all_loc={};ID=[];
    switch batchi
        case 1
            fishind=numL;
        case 4
            fishind=numC;
        case 3
            fishind=numNL;
        case 2
            fishind=numFL;
    end
    for typei=1:10
        f=['ind_all_CSUS_RESPONSIVE',labels_typei_ADJ{typei}];a=getfield(A,f);
        for sessioni=1:size(a,1)
            for ii=1:length(fishind)
                p=path_all{fishind(ii)};
                nn=[p(end-14:end-7),p(end-5:end-1)];
                supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                ind=a(sessioni,fishind(ii));
                loc{ii}=supervoxeli(id,1:3);
                all_loc{ii}=supervolxeli(:,1:3);
                ID(ii,1:length(id))=id;
            end
            [P_z,H_z,P_t,H_t,Dr_ij,Dn_ij]=statistic_P_in_spatical_loc(loc,all_loc);
            save([checkpath(fullfile(savepath,'p_spatial-types\')),[seg_batchi{batchi},labels_typei_ADJ{typei},sessionx(sessioni)],'-P_spatial.mat'],'loc','all_loc','ID', 'Dn_ij','Dr_ij', 'P_z','H_z','P_t','H_t', '-v7.3');
        end
    end
end
