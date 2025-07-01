clc;clear all;
res=[0.66,0.66,10];radium=floor(30);

load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
%id=find(temp_supervoxel(:,3)==(30-1)*res(3));figure,scatter(temp_supervoxel(id,1),temp_supervoxel(id,2))
load('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\Path.mat');
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
mat_name=[];
seg_batchi={'Learner','Unpair','Non-learner','Faded-learner'}; %seg_batchi={'learner','Unpair','Non-learner','Faded-learner'}; 
seg1={'Activation','Inhibition','Unresponsive'};%seg1={'ON','OFF','UNRESPONSIVE'};
seg2={'Increased','Decreased','Emerged','Shift to other','Increased2'};clr=jet(length(seg1)*length(seg2));
Fraction_type={};
%% pre_post
%学习前后CR pattern
CR_ind_all_pre=cell(4,2);CR_ind_all_post=cell(4,2);
fraction={};Fraction_in_region_pre={};Fraction_in_region_post={};
for batchi=1:length(Path)
    path=Path{batchi};
    for ii=1:length(path)
        load([path{ii} 'CR_ind_summary.mat']);
        nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)]
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        savepath_CR_pattern=checkpath(fullfile(savepath,'CR pattern'));
        if ~isempty(strfind(nn,'2019'))
            res2=[0.66,0.66,8];
        elseif ~isempty(strfind(nn,'2021')) | ~isempty(strfind(nn,'2022'))
            res2=[0.66,0.66,10];
        end
        for session=1:2
            switch session
                case 1 %pre
                    for type=1:2
                        switch type
                            case 1
                                ind=find(CR_ind_up(:,1)==1);
                                %ind=intersect(ind,CS_related_hab_ind);
                            case 2
                                ind=find(CR_ind_down(:,1)==1);
                                %ind=intersect(ind,CS_related_hab_ind);
                        end
                    %temp_plot( supervolxeli(ind,1:3),temp_supervoxel,temp_env.env.vol,['pre Cond. CR_',nn],1, savepath_CR_pattern,nn,[])
                    fraction{batchi,type}(ii,1)=length(ind)/size(CR_ind_up,1);
                    [In,~]=count_spatial_location_density_temp(supervolxeli(ind,1:3),temp_supervoxel, radium,nn);
                    %temp_plot( supervolxeli(ind,1:3),temp_supervoxel,temp_env.env.vol,['pre Cond. CR-',nn,'-clrmap'],2, savepath_CR_pattern,nn,In)
                    CR_ind=[supervolxeli(ind,:),string(repmat(nn,length(ind),1)),In];
                    CR_ind_all_pre{batchi,type}=[CR_ind_all_pre{batchi,type};CR_ind];
                    
                    surpervolxel=supervolxeli(ind,1:3);%double([CR_ind_all_pre{batchi}(:,1:3);CR_ind_all_post{batchi}(:,1:3)]);
                    gIX=[ones(size( surpervolxel,1),1)];%[ones(length(CR_ind_all_pre{batchi}),1);2*ones(length(CR_ind_all_post{batchi}),1)];
                    gIX(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];surpervolxel(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
                    gIX(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];surpervolxel(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
                    gIX(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)))=[];surpervolxel(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)),:)=[];
                    [~,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
                        gIX,surpervolxel,temp_supervoxel,clr,false,nn);
                    Fraction_in_region_pre{batchi,type}(:,ii)=fraction_in_region_in_clust(:,1);
                    end
                case 2  %post
                    for type=1:2
                        switch type
                            case 1
                                ind=find(CR_ind_up(:,2)==1);
                                %ind=intersect(ind,CS_related_tst_ind);
                            case 2
                                ind=find(CR_ind_down(:,2)==1);
                                %ind=intersect(ind,CS_related_tst_ind);
                        end
                        %temp_plot( supervolxeli(ind,:),temp_supervoxel,temp_env.env.vol,['post Cond. CR_',nn],1, savepath_CR_pattern,nn,[])
                        fraction{batchi,type}(ii,2)=length(ind)/size(CR_ind_up,1);
                        [In,~]=count_spatial_location_density_temp(supervolxeli(ind,1:3),temp_supervoxel, radium,nn);
                        %temp_plot( supervolxeli(ind,1:3),temp_supervoxel,temp_env.env.vol,['post Cond. CR-',nn,'-clrmap'],2, savepath_CR_pattern,nn,In)
                        CR_ind=[supervolxeli(ind,:),string(repmat(nn,length(ind),1)),In];
                        CR_ind_all_post{batchi,type}=[CR_ind_all_post{batchi,type};CR_ind];
                        
                        surpervolxel=supervolxeli(ind,1:3);%double([CR_ind_all_pre{batchi}(:,1:3);CR_ind_all_post{batchi}(:,1:3)]);
                        gIX=[ones(size( surpervolxel,1),1)];%[ones(length(CR_ind_all_pre{batchi}),1);2*ones(length(CR_ind_all_post{batchi}),1)];
                        gIX(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];surpervolxel(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
                        gIX(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];surpervolxel(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
                        gIX(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>size(temp_env.env.vol,3)*res2(3)))=[];surpervolxel(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>size(temp_env.env.vol,3)*res2(3)),:)=[];
                        [~,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
                            gIX,surpervolxel,temp_supervoxel,clr,false,nn);
                        Fraction_in_region_post{batchi,type}(:,ii)=fraction_in_region_in_clust(:,1);
                    end
            end
        end
    end
end
save([savepath,'fraction_CR_pre_post.mat'],'CR_ind_all_pre','CR_ind_all_post','fraction','Fraction_in_region_pre','Fraction_in_region_post', '-v7.3');
%% CR_stastic.mat
for batchi=1:length(Path)
    path=Path{batchi};
    for ii=1:length(path)
        load([path{ii} 'CR_ind_summary.mat']);
        nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)];
        A=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        if 1%~exist(fullfile(path{ii},'CR_stastic.mat'))
            %load(fullfile(path{ii},'CR_stastic.mat'));
            m_pre=squeeze(mean(CS_response{1},1,'omitnan'));m_aft=squeeze(mean(CS_response{2},1,'omitnan'));
            p_increase=nan(size(m_pre,2),1);p_decrease=nan(size(m_pre,2),1);
            p_increase_nonabs=nan(size(m_pre,2),1);p_decrease_nonabs=nan(size(m_pre,2),1);
            parfor jj=1:size(m_pre,2)
                a=m_pre(:,jj);a(isnan(a))=[];b=m_aft(:,jj);b(isnan(b))=[];
                p_increase(jj)=ranksum(abs(a),abs(b),'tail','left');%a<b
                p_decrease(jj)=ranksum(abs(a),abs(b),'tail','right');%p<0.05,接受H1:A>B
                p_increase_nonabs(jj)=ranksum(a,b,'tail','left');%a<b
                p_decrease_nonabs(jj)=ranksum(a,b,'tail','right');%p<0.05,接受H1:A>B
            end
            save(fullfile(path{ii},'CR_stastic.mat'),'p_increase','p_decrease','p_increase_nonabs','p_decrease_nonabs', '-v7.3');
        else
            load(fullfile(path{ii},'CR_stastic.mat'));
        end
    end
end

%% types_pre_post
%load(fullfile(savepath,'fraction_types.mat'));
Fraction_in_region_type={}; Fraction_type={}; Fraction_in_region_type_in_cluster={};
for batchi=1:length(Path)
    path=Path{batchi};
    for zz=1:length(seg1)
        for hh=1:length(seg2)
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            if exist(csv_file)
                delete(csv_file);
            end
            MM=table;
            for ii=1:length(path)
                load([path{ii} 'CR_ind_summary.mat']);
                nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)];
                if ~isempty(strfind(nn,'2019'))
                    res2=[0.66,0.66,8];
                elseif ~isempty(strfind(nn,'2021'))
                    res2=[0.66,0.66,10];
                end
                A=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                load(fullfile(path{ii},'CR_stastic.mat'));
                switch seg1{zz}
                    case 'ON'
                        ind=find(CR_ind_up(:,1)==1);
                    case 'OFF'
                        ind=find(CR_ind_down(:,1)==1);
                    case 'UNRESPONSIVE'
                        ind=find(CR_ind_up(:,1)==1);
                        %ind=intersect(ind,CS_related_hab_ind);
                        ind2=find(CR_ind_down(:,1)==1);
                        %ind2=intersect(ind2,CS_related_hab_down_ind);
                        ind3=union(ind,ind2);
                        ind=setdiff(1:size(CR_ind_up(:,1),1),ind3);
                end
                switch seg2{hh}
                    case 'Increased'
                        if zz==3
                            ind2=find(p_increase_nonabs<0.05);
                            ind3=find(CR_ind_up(:,2)==1);
                            a=intersect(ind2,ind3);
                            IND=intersect(a,ind);
                        else
                            ind2=find(p_increase<0.05);
                            IND=intersect(ind2,ind);
                        end
                    case 'Increased2'
                        ind2=find(p_increase<0.05);
                        IND=intersect(ind2,ind);
                    case 'Decreased'
                        if zz==3
                            ind2=find(p_decrease_nonabs<0.05);
                            ind3=find(CR_ind_down(:,2)==1);
                            a=intersect(ind2,ind3);
                            IND=intersect(a,ind);
                        else
                        ind2=find(p_decrease<0.05);
                        IND=intersect(ind2,ind);
                        end
                    case 'Emerged'
                        if zz==1 
                            ind=find(CR_ind_up(:,1)==1);
                            ind2=find(CR_ind_down(:,1)==1);
                            ind3=union(ind,ind2);
                            ind=setdiff(1:size(CR_ind_up(:,1),1),ind3);
                            ind2=find(CR_ind_up(:,2)==1);
                            ind3=find(p_increase<0.05);
                            a=intersect(ind2,ind3);
                            IND=intersect(a,ind);
                        elseif zz==2
                            ind=find(CR_ind_up(:,1)==1);
                            ind2=find(CR_ind_down(:,1)==1);
                            ind3=union(ind,ind2);
                            ind=setdiff(1:size(CR_ind_up(:,1),1),ind3);
                            ind2=find(CR_ind_down(:,2)==1);
                            ind3=find(p_increase<0.05);
                            a=intersect(ind2,ind3);
                            IND=intersect(a,ind);
                        end
                        %ind2=intersect(ind2,CS_related_tst_ind);
                        %IND=intersect(ind2,ind);
                    case 'shift to other'
                        if zz==1
                            ind2=find(CR_ind_down(:,2)==1);
                            %ind2=intersect(ind2, CS_related_tst_down_ind);
                            IND=intersect(ind2,ind);
                        elseif zz==2
                            ind2=find(CR_ind_up(:,2)==1);
                            %ind2=intersect(ind2,CS_related_tst_ind);
                            IND=intersect(ind2,ind);
                        elseif zz==3
                            IND=[];
                        end
                end
                x=A(IND,1);y=A(IND,2);z=A(IND,3);
                columns = {'x', 'y', 'z','ind','label'};
                if 0<length(IND) & length(IND)<=1
                    M=table(x,y,z,IND,categorical(cellstr(repmat(nn,length(IND),1))),'VariableNames', columns);
                elseif length(IND)>1
                    M=table(x,y,z,IND,categorical(cellstr(repmat(nn,length(IND),1))),'VariableNames', columns);
                    %M=table(x,y,z,IND,mat2cell(repmat(nn,length(IND),1),ones(length(IND),1),length(nn)*ones(length(IND),1)),'VariableNames', columns);
                end
                MM=cat(1,MM,M);%[MM;M];
                Fraction_type{batchi,zz}(hh,ii)=length(IND)./length(CR_ind_up(:,1));
                
                surpervolxel=[x,y,z];
                gIX=[ones(size( surpervolxel,1),1)];
                gIX(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];surpervolxel(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
                gIX(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];surpervolxel(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
                gIX(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)))=[];surpervolxel(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)),:)=[];
                [~,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
                    gIX,surpervolxel,temp_supervoxel,clr,false,nn);
                Fraction_in_region_type{batchi,zz}(:,hh,ii)=fraction_in_region_in_clust(:,1);
                Fraction_in_region_type_in_cluster{1}{batchi,zz}(:,hh,ii)= num_in_region_in_clust{1}(:,1);
                Fraction_in_region_type_in_cluster{2}{batchi,zz}(:,hh,ii)= num_in_region_in_clust{2}(:,1);
                Fraction_in_region_type_in_cluster{3}{batchi,zz}(:,hh,ii)= num_in_region_in_clust{3}(:,1);
            end
            writetable(MM,csv_file,'WriteVariableNames',true,'WriteRowNames',false);
            %writetable(M,csv_file,'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);
        end
    end
end
save([savepath,'fraction_types.mat'],'Fraction_in_region_type','Fraction_type','Fraction_in_region_type_in_cluster', '-v7.3');

%% all_loc
for batchi=1:length(Path)
    all_loc=table;
    path=Path{batchi};
    for ii=1:length(path)
        nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)]
        supervolxeli=readmatrix([path{ii} ,nn,'vol_cut_env_spatialloc_warped_SyN.csv']);
        x=supervolxeli(:,1);y=supervolxeli(:,2);z=supervolxeli(:,3);IND=supervolxeli(:,5);
        columns = {'x', 'y', 'z','ind','label'};
        if 0<length(IND) & length(IND)<=1
            M=table(x,y,z,IND,categorical(cellstr(repmat(nn,length(IND),1))),'VariableNames', columns);
        elseif length(IND)>1
            M=table(x,y,z,IND,categorical(cellstr(repmat(nn,length(IND),1))),'VariableNames', columns);
            %M=table(x,y,z,IND,mat2cell(repmat(nn,length(IND),1),ones(length(IND),1),length(nn)*ones(length(IND),1)),'VariableNames', columns);
        end
        all_loc=cat(1,all_loc,M);%[MM;M];
    end
    save(fullfile(savepath,[seg_batchi{batchi},'-all_loc.mat']),'all_loc', '-v7.3');
end

%% 统计X-Y分布  
load(fullfile(savepath,'fraction_CR_pre_post.mat'))
%CR_pre_post
nbin=size(temp_env.env.vol);
X_Y_Z_Density={};
for batchi=1:length(Path)
    for CR_type=1:2
        loc=double(CR_ind_all_pre{batchi,CR_type}(:,1:3)); %figure;
        for type=1:3 %x,y,z
            a=loc(:,type);
            h=histogram(a,nbin(type),'visible','off');
            x=h.BinEdges(1:end-1)'; y=h.Values';z=ones(length(x),1)*length(a);
            X_Y_Z_Density{batchi,CR_type,type}=[x,y,z];
            %subplot(3,1,type),plot(x,y,'k','linewidth',1.5);box off;set(gca,'fontsize',14);
        end
    end
end
save([savepath,'density_CR_pre_post.mat'],'X_Y_Z_Density', '-v7.3');
%CR_types
X_Y_Z_Density={};
for batchi=1:length(Path)
    for zz=1:3
        for hh=[1 2 4]
            if strcmp(seg1{zz},'UNRESPONSIVE') && strcmp(seg2{hh},'shift to other')
                continue;
            end
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);loc= [M.x M.y M.z];
            for type=1:3 %x,y,z
                a=loc(:,type);
                h=histogram(a,nbin(type),'visible','off');
                x=h.BinEdges(1:end-1)'; y=h.Values';z=ones(length(x),1)*length(a);
                X_Y_Z_Density{batchi,zz,hh,type}=[x,y,z];
            end
        end
    end
end
save([savepath,'density_CR_type.mat'],'X_Y_Z_Density', '-v7.3');
%
%% 手动统计
Lable={'L OB','R OB','L P','R P',...
    'L Hb','R Hb','L TeO','R TeO','L Np','R Np','L Otr','R Otr','TL','PO',...
    'Pr','PT','Th','rH','iH','cH','mPT','T',...
    'L TS','R TS','Va','CCe','L gT','R gT','R1','R2'};
% for batchi=1:4
%     path=Path{batchi};
%     for zz=1:length(seg1)
%         for hh=1:length(seg2)
%             for ii=1:length(path)
%
%             end
%         end
%     end
% end
load('H:\3.Juvenile reference brain\registration to templete\data\2021_ind\fraction_CR_pre_post.mat')
type=2;
A=[];
for zz=1:2
    aa=[];
    for batchi=1:3
        a=(fraction{batchi,type}(:,zz))';
        if length(a)<6
            a=[a,nan,nan];
        end
        aa=[aa,a];
    end
    A=[A;aa];
end
p=[];
for batchi=1:3
    a=(fraction{batchi,type}(:,1))';
    b=(fraction{batchi,type}(:,2))';
    p(batchi,type)=ranksum(a,b);
end
%% Fraction_type
load('H:\3.Juvenile reference brain\registration to templete\data\2021_ind\fraction_types.mat')
kk=1;
A=[];
for zz=1:3
    for hh=[1 2 4]
        aa=[];
        for batchi=1:3
            a=(Fraction_type{batchi,zz}(hh,:));
            if length(a)<6
                a=[a,nan,nan];
            end
            aa=[aa,a];
        end
        A(kk,:)=aa;kk=kk+1;
    end
end
p=[];
for ii=1:3
    for jj=1:3
        for zz=1:3
            for hh=[1 2 3 4]             
                a=(Fraction_type{ii,zz}(hh,:));
                b=(Fraction_type{jj,zz}(hh,:));
                p(ii,jj,zz,hh)=ranksum(a,b);
            end
        end
    end
end
[x,y,z,t]=ind2sub(size(p),find(p<0.05));
for ii=1:length(x)
    disp([seg_batchi{x(ii)},'-',seg_batchi{y(ii)},seg1{z(ii)},seg2{t(ii)}])
end
%% Fraction_type_region  
%脑区各类占比 /脑区内总neuron
aa=[];b=nan(30,2);
for zz=1:3
    for hh=[1 2 4]
        a=(squeeze(Fraction_in_region_type{3,zz}(:,hh,:))*100);
        if size(a,2)<6
            a=[a,b];
        end
       aa=[aa,a];
    end
    
end
%脑区各类占比-转置格式
b=nan(1,2);aaa=[];
for zz=1:3
    for hh=[1 2 4]
        aa=[];
        for regioni=1:30
        a=(squeeze(Fraction_in_region_type{1,zz}(regioni,hh,:))*100)';
        if length(a)<6
            a=[a,b];
        end
        aa=[aa,a];
        end
        aaa=[aaa;aa];
    end
end
%脑区内各类占比 /脑区内engram cell neuron
bb=[];
for ii=1:6
    aa=[];
    for zz=1:3
        for hh=[1 2 4]
            kk=Fraction_in_region_type{3,zz};
            if ii>size(kk,3)
                a=nan(30,1);
            else
               a=(squeeze(kk(:,hh,ii))*100);
            end
            aa=[aa,a];
        end
    end
    aaa=aa./repmat(sum(aa,2,'omitnan'),1,9);
    bb(:,ii:6:6*9)=aaa;
end
%每类在脑区的占比
aa=[];b=nan(30,2);
for ii=1:3
    a=squeeze(Fraction_in_region_type{ii,3}(:,2,:))*100;
    if size(a,2)<6
        a=[a,b];
    end
    aa=[aa,a];
end
p=ones(3,3,3,4,30);
for ii=1:3
    for jj=setdiff([1:3],ii)
        for zz=1:3
            for hh=[1 2 4]
                for region=1:30
                a=squeeze(Fraction_in_region_type{ii,zz}(region,hh,:));
                b=squeeze(Fraction_in_region_type{jj,zz}(region,hh,:));
                if sum(isnan(a))==length(a)
                  p(ii,jj,zz,hh,region)=nan;
                else
                    p(ii,jj,zz,hh,region)=ranksum(a,b);
                end
                end
            end
        end
    end
end
[x,y,z,t,h]=ind2sub(size(p),find(p<0.05));
for ii=1:length(x)
    disp([seg_batchi{x(ii)},'-',seg_batchi{y(ii)},seg1{z(ii)},seg2{t(ii)},'-',Lable{h(ii)}])
end
%% Fraction_type_P_spatial_0.05
load([checkpath(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\')),'Fraction_in_region_type_p_spatial.mat']);
 kk=1;
A=[];
for zz=1:3
    for hh=[1 2 4]
        aa=[];
        if zz==3 & hh==4
            continue;
        end
        for batchi=1:3
            a=(Fraction_type_p_spatial{batchi,zz}(hh));
            aa=[aa,a];
        end
        A(kk,:)=aa;kk=kk+1;
    end
end
%% Fraction_type_region_P_spatial_0.05
%脑区各类占比 /脑区内总neuron
aa=[];b=nan(30,2);
for zz=1:3
    for hh=[1 2 4]
        if zz==3 && hh==4
            continue;
        end
        a=(squeeze(Fraction_in_region_type_p_spatial_fishi{3,zz,1}(:,hh,:))*100);
        if size(a,2)<6
            a=[a,b];
        end
        aa=[aa,a];
    end
end
%脑区内各类占比 /脑区内engram cell neuron
bb=[];
for ii=1:6
    aa=[];
    for zz=1:3
        for hh=[1 2 4]
            if zz==3 && hh==4
                continue;
            end
            kk=Fraction_in_region_type_p_spatial_fishi{3,zz,1};
            if ii>size(kk,3)
                a=nan(30,1);
            else
               a=(squeeze(kk(:,hh,ii))*100);
            end
            aa=[aa,a];
        end
    end
    aaa=aa./repmat(sum(aa,2,'omitnan'),1,8);
    bb(:,ii:6:6*8)=aaa;
end
%每类在脑区的占比
aa=[];b=nan(30,2);
for ii=1:3
    a=squeeze(Fraction_in_region_type_p_spatial_fishi{ii,3,1}(:,1,:))*100;
    if size(a,2)<6
        a=[a,b];
    end
    aa=[aa,a];
end
p=ones(3,3,3,4,30);
for ii=1:3
    for jj=setdiff([1:3],ii)
        for zz=1:3
            for hh=[1 2 4]
                if zz==3 && hh==4
                    continue;
                end
                for region=1:30
                a=squeeze(Fraction_in_region_type_p_spatial_fishi{ii,zz,1}(region,hh,:));
                b=squeeze(Fraction_in_region_type_p_spatial_fishi{jj,zz,1}(region,hh,:));
                if sum(isnan(a))==length(a)
                  p(ii,jj,zz,hh,region)=nan;
                else
                    p(ii,jj,zz,hh,region)=ranksum(a,b);
                end
                end
            end
        end
    end
end
[x,y,z,t,h]=ind2sub(size(p),find(p<0.05));
for ii=1:length(x)
    disp([seg_batchi{x(ii)},'-',seg_batchi{y(ii)},seg1{z(ii)},seg2{t(ii)},'-',Lable{h(ii)}])
end
%all 
aa=[];
for ii=1:3
    a=squeeze(Fraction_in_region_type_p_spatial{ii,3,1}(:,2))*100;
    aa=[aa,a];
end
%% Fraction_in_region_pre_post
type=2;
aa=[];
for ii=1:3
    a=Fraction_in_region_pre{ii,type};
    b=Fraction_in_region_post{ii,type};
    aa=[a,b];
    if size(a,2)<6
        a=[a,nan(size(a,1),2)];
    end
    aa=[aa,a];
end
p=[];
for jj=1:3
    for zz=1:3
        for ii=1:30
            a=Fraction_in_region_post{jj,type}(ii,:);
            b=Fraction_in_region_post{zz,type}(ii,:);
            p(jj,zz,ii)=ranksum(a,b);
        end
    end
end
[x,y,z]=ind2sub(size(p),find(p<0.05));
for ii=1:length(x)
    disp([seg_batchi{x(ii)},'-',seg_batchi{y(ii)},' ',Lable{z(ii)}])
end

for zz=1:3
    p=[];
    for ii=1:30
        a=Fraction_in_region_pre{zz,type}(ii,:);
        b=Fraction_in_region_post{zz,type}(ii,:);
        p(ii)=ranksum(a,b);
        %[~,p(ii)] = ttest2(a,b);
    end
    x=find(p<0.05);
    for ii=1:length(x)
        disp([seg_batchi{zz},Lable{x(ii)}])
    end
end

%% re-organize figure
figure('position',[7,154,1910,824]),t=tiledlayout(8,3);t.TileSpacing = 'compact';t.Padding = 'compact';
zz=[1,1;3,1;1,2;1,4;2,1;3,2;2,2;2,4];
for jj=1:size(zz,1)
    for batchi=1:3
        pp=fullfile(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\'),['Max_',[seg_batchi{batchi},seg1{zz(jj,1)},seg2{zz(jj,2)},'-P_spatial','-type',num2str(1)],'.tif']);
        I=imread(pp);
        nexttile,imshow(I)
    end
end
