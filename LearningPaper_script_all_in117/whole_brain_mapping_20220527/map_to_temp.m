
savepath='H:\3.Juvenile reference brain\registration to templete\data\2021_ind\';
%% trace of differnt type
for batchi=1
    path=Path{batchi};aa=string(path)';aa=strrep(aa,'\','');
    for zz=1:length(seg1)
        for hh=1:length(seg2)
            CR_act_pre={};CR_act_post={};
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);supervoxel= [M.x M.y M.z];a=length(string(unique(M.label))');kk=1; figure('name',[seg_batchi{batchi},seg1{zz},seg2{hh}]),         
            for name=string(unique(M.label))'
               fish_i= find(contains(aa,name));p=path{fish_i};
               ind_in_fishi=M.ind(find(contains(M.label,name)));
               load(fullfile(p,'activities_aft_process.mat'));load(fullfile(p,'para.mat'));load(fullfile(p,'env.mat'));
               trial_hab=[ trial.hab(2):trial.hab(3)];trial_test=[trial.test(2):trial.test(3)];
               A=activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,ind_in_fishi);A_r=reshape(A,frame.per_cycle,trial.total,[]);
               aa_pre=squeeze(mean(A_r(:,trial_hab,:),2));aa_post=squeeze(mean(A_r(:,trial_test,:),2));
               subplot(2,ceil(a/2),kk);plot(mean(aa_pre,2),'k--','linewidth',2);hold on;plot(mean(aa_post,2),'k','linewidth',2);ylabel('ΔF/F');set(gca,'fontsize',16);title(name);
               kk=kk+1;
               CR_act_pre{fish_i}=aa_pre;CR_act_post{fish_i}=aa_post;
            end    
            save([checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\','CR_act_pre_post','.mat'],'CR_act_pre','CR_act_post', '-v7.3');
        end
    end
end
for batchi=1:3
    path=Path{batchi};aa=string(path)';aa=strrep(aa,'\','');
    for zz=1:length(seg1)
        for hh=1:length(seg2)
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);supervoxel= [M.x M.y M.z];a=length(string(unique(M.label))');kk=1; figure('name',[seg_batchi{batchi},seg1{zz},seg2{hh}]),
            load([checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\','CR_act_pre_post','.mat']);
            for name=string(unique(M.label))'
               fish_i= find(contains(aa,name));p=path{fish_i};
               aa_pre = CR_act_pre{fish_i};aa_post=CR_act_post{fish_i};
               subplot(2,ceil(a/2),kk);plot(mean(aa_pre,2),'k--','linewidth',2);hold on;plot(mean(aa_post,2),'k','linewidth',2);ylabel('ΔF/F');set(gca,'fontsize',16);title(name);
               kk=kk+1;
            end       
        end
    end
end
%% map different type
is_plot=true;
count_spatial_p=true;
count_fraction_p=false;
for batchi=1:3
    for zz=[3]
        for hh=[1 2]
            if strcmp(seg1{zz},'UNRESPONSIVE') && strcmp(seg2{hh},'shift to other')
                continue;
            end
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);supervoxel= [M.x M.y M.z];
            %plot
            if is_plot
            temp_plot( supervoxel,temp_env.env.supervoxel(:,1:3),temp_env.env.vol,[seg_batchi{batchi},seg1{zz},seg2{hh}],1,checkpath(fullfile(savepath,'p_spatial-types\')),'2021',[]);
            end
            %p value in_spatical
            if count_spatial_p
                path=Path{batchi};loc={};all_loc={};ID=[];
                for ii=1:length(path)
                    nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)]
                    id=find(strcmp( M.label,nn));
                    loc{ii}=supervoxel(id,1:3);
                    supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                    all_loc{ii}=supervolxeli(:,1:3);
                    ID(ii,1:length(id))=id;
                end
                [P_z,H_z,P_t,H_t,Dr_ij,Dn_ij]=statistic_P_in_spatical_loc(loc,all_loc);
                save([checkpath(fullfile(savepath,'p_spatial-types\')),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_spatial.mat'],'loc','all_loc','ID', 'Dn_ij','Dr_ij', 'P_z','H_z','P_t','H_t', '-v7.3');
            end
            %p value in_fraction of region
            if count_fraction_p
                load(fullfile(savepath,[seg_batchi{batchi},'-all_loc.mat']));
                load(fullfile(savepath,'fraction_types.mat'));
                [P_z,H_z,P_t,H_t,fraction_samlpe]=statistic_P_in_region_fraction(M,squeeze(Fraction_in_region_type{batchi,zz}(:,hh,:)),all_loc,temp_env,temp_supervoxel,reg_mask,reg_name,reg_loc);
                save([checkpath(fullfile(savepath,'p_fraction-types\')),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_fraction.mat'],'all_loc','Fraction_in_region_type','fraction_samlpe', 'P_z','H_z','P_t','H_t', '-v7.3');
            end
        end
    end
end
%% 重新计算p_spatial-types
for batchi=1:3
    for zz=[3]
        for hh=[1 2]
            pp=[fullfile(savepath,'p_spatial-types\'),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_spatial.mat']
            if exist(pp)
            load(pp);
            [P_z,H_z,Dr_ij_projection,Dn_ij_projection]=adjust_P_spatial_loc(loc,Dr_ij,Dn_ij);
            save([checkpath(fullfile(savepath,'p_spatial-types\')),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_spatial_adjust_projectioned_60.mat'],...
                'loc','all_loc','ID', 'Dn_ij','Dr_ij', 'P_z','H_z','P_t','H_t', 'Dr_ij_projection','Dn_ij_projection','-v7.3');
            end
        end
    end
end
%% plot p value-types
Fraction_in_region_type_p_spatial={};Fraction_in_region_type_in_cluster_p_spatial={};
Fraction_in_region_type_p_spatial_fishi={};Fraction_in_region_type_in_cluster_p_spatial_fishi={};
Fraction_type_p_spatial={};Fraction_type_p_spatial_fishi={};
for batchi=1:3
    for zz=1:length(seg1)
        for hh=1:length(seg2)
            ppath=[fullfile(savepath,'p_spatial-types\'),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_spatial_adjust_projectioned_60.mat'];
            if exist(ppath)
                csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
                M=readtable(csv_file);supervoxel= [M.x M.y M.z];label=M.label;
                load(ppath);
                P=P_z;P_z_r=reshape(P,1,[]); ID_r=reshape(ID,1,[]);
                P_z_r(find(ID_r==0))=[];ID_r(find(ID_r==0))=[];
                for type=1
                    switch type
                        case 1
                            In=(1-P_z_r);P_r=P_z_r;
                    end
                    kkk=1;
                    for pp=[0.05,0.025]
                        ind=find(P_r<pp);
                        temp_plot( supervoxel(ind,:),temp_env.env.supervoxel(:,1:3),temp_env.env.vol,[seg_batchi{batchi},seg1{zz},seg2{hh},'-P_spatial_',num2str(pp),'-type',num2str(type)],1,checkpath(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\')),'2021',[])
                        surpervolxel_p=supervoxel(ind,:);fish=label(ind);kkkk=1;
                        Fraction_type_p_spatial{batchi,zz}(hh)=length(ind)./size(temp_supervoxel,1);
                        gIX=[ones(size( surpervolxel_p,1),1)];res2=[0.66,0.66,10];
                        gIX(find(surpervolxel_p(:,1)<0 | surpervolxel_p(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];surpervolxel_p(find(surpervolxel_p(:,1)<0 | surpervolxel_p(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
                        gIX(find(surpervolxel_p(:,2)<0 | surpervolxel_p(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];surpervolxel_p(find(surpervolxel_p(:,2)<0 | surpervolxel_p(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
                        gIX(find(surpervolxel_p(:,3)<0 | surpervolxel_p(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)))=[];surpervolxel_p(find(surpervolxel_p(:,3)<0 | surpervolxel_p(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)),:)=[];
                        [~,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
                            gIX,surpervolxel_p,temp_supervoxel,clr,false,'2021');
                        Fraction_in_region_type_p_spatial{batchi,zz,kkk}(:,hh)=fraction_in_region_in_clust(:,1);
                        Fraction_in_region_type_in_cluster_p_spatial{1}{batchi,zz,kkk}(:,hh)= num_in_region_in_clust{1}(:,1);
                         Fraction_in_region_type_in_cluster_p_spatial{2}{batchi,zz,kkk}(:,hh)= num_in_region_in_clust{2}(:,1);
                          Fraction_in_region_type_in_cluster_p_spatial{3}{batchi,zz,kkk}(:,hh)= num_in_region_in_clust{3}(:,1);
                        surpervolxel_p=supervoxel(ind,:);
                        for ii=unique(fish)'
                            if ~isempty(strfind(ii,'2019'))
                                res2=[0.66,0.66,8];
                            elseif ~isempty(strfind(ii,'2021'))
                                res2=[0.66,0.66,10];
                            end
                            pppp=find(contains(strrep(string(Path{batchi}),'\',''),ii));
                            load(string(fullfile(Path{batchi}(pppp),'env.mat')));
                            fish_i= find(strcmp(ii,fish));surpervolxel=surpervolxel_p(fish_i,:);gIX=[ones(size( surpervolxel,1),1)];
                            Fraction_type_p_spatial_fishi{batchi,zz}(hh,kkkk)=length(fish_i)./size(env.supervoxel,1);
                            gIX(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];surpervolxel(find(surpervolxel(:,1)<0 | surpervolxel(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
                            gIX(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];surpervolxel(find(surpervolxel(:,2)<0 | surpervolxel(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
                            gIX(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)))=[];surpervolxel(find(surpervolxel(:,3)<0 | surpervolxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)),:)=[];
                            if ~isempty(surpervolxel)
                            [~,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
                                gIX,surpervolxel,temp_supervoxel,clr,false,ii);
                            Fraction_in_region_type_p_spatial_fishi{batchi,zz,kkk}(:,hh,kkkk)=fraction_in_region_in_clust(:,1);
                            Fraction_in_region_type_in_cluster_p_spatial_fishi{1}{batchi,zz,kkk}(:,hh,kkkk)= num_in_region_in_clust{1}(:,1);
                            Fraction_in_region_type_in_cluster_p_spatial_fishi{2}{batchi,zz,kkk}(:,hh,kkkk)= num_in_region_in_clust{2}(:,1);
                            Fraction_in_region_type_in_cluster_p_spatial_fishi{3}{batchi,zz,kkk}(:,hh,kkkk)= num_in_region_in_clust{3}(:,1);
                            kkkk=kkkk+1; 
                            end
                        end
                    end
                    kkk=kkk+1;
                    In=floor(In*100);
                    temp_plot( supervoxel,temp_env.env.supervoxel(:,1:3),temp_env.env.vol,[seg_batchi{batchi},seg1{zz},seg2{hh},'-P_spatial','-type',num2str(type)],2,checkpath(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\')),'2021',In);
                end
            end
        end
    end
end
save([checkpath(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\')),'Fraction_in_region_type_p_spatial.mat'],'Fraction_type_p_spatial','Fraction_type_p_spatial_fishi',...
    'Fraction_in_region_type_p_spatial','Fraction_in_region_type_in_cluster_p_spatial','Fraction_in_region_type_p_spatial_fishi','Fraction_in_region_type_in_cluster_p_spatial_fishi','-v7.3');
%% statistic of Fraction_type_p_spatial
A=Fraction_in_region_type_p_spatial;
for batchi=1:size(A,1)
    for zz=1:size(A,2)
        for hh=1:size(A{batchi,zz},2)
            csv_file=[checkpath(fullfile(savepath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);supervoxel= [M.x M.y M.z];

            all_loc=table(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),[1:size(temp_supervoxel,1)]',categorical(cellstr(repmat('2021',size(temp_supervoxel,1),1))),'VariableNames', {'x', 'y', 'z','ind','label'});
            [P_z,H_z,P_t,H_t,fraction_samlpe]=statistic_P_in_region_fraction(M,squeeze(Fraction_in_region_type_p_spatial{batchi,zz,1}(:,hh)),all_loc,temp_env,temp_supervoxel,reg_mask,reg_name,reg_loc);
            save([checkpath(fullfile(savepath,'p_spatial-types_adjust_projectioned_60\')),[seg_batchi{batchi},seg1{zz},seg2{hh}],'-P_fraction_p_spatial.mat'],'all_loc','Fraction_in_region_type','fraction_samlpe', 'P_z','H_z','P_t','H_t', '-v7.3');
        end
    end
end
%% plot pre_post CR
load([savepath,'fraction_CR_pre_post.mat']);
for type=1:2
    for batchi=1:3
        supervoxel= double(CR_ind_all_pre{batchi,type}(:,1:3));
        In=double(CR_ind_all_pre{batchi,type}(:,7));
        temp_plot( supervoxel,temp_env.env.supervoxel(:,1:3),temp_env.env.vol,['pre Cond. CR-',seg_batchi{batchi},'-type',num2str(type)],2,savepath,'2021',In)
        
        supervoxel= double(CR_ind_all_post{batchi,type}(:,1:3));
        In=double(CR_ind_all_post{batchi,type}(:,7));
        temp_plot( supervoxel,temp_env.env.supervoxel(:,1:3),temp_env.env.vol,['post Cond. CR-',seg_batchi{batchi},'-type',num2str(type)],2,savepath,'2021',In)
    end
end
%% 计算pre_post CR spatial_loc
for type=1:2
    for batchi=1
        path=Path{batchi};
        loc_pre={};all_loc={};loc_post={};
        ID_pre=[];ID_post=[];
        %
        for ii=1:length(path)
            nn=[path{ii}(end-14:end-7),path{ii}(end-5:end-1)]
            supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
            id_pre=find(strcmp( CR_ind_all_pre{batchi,type}(:,6),nn));
            loc_pre{ii}=double(CR_ind_all_pre{batchi,type}(id_pre,1:3));
            
            id_post=find(strcmp( CR_ind_all_post{batchi,type}(:,6),nn));
            loc_post{ii}=double(CR_ind_all_post{batchi,type}(id_post,1:3));
            all_loc{ii}=supervolxeli(:,1:3);
            ID_pre(ii,1:length(id_pre))=id_pre;
            ID_post(ii,1:length(id_post))=id_post;
        end
        
        [P_pre_z,H_pre_z,P_pre_t,H_pre_t,Dr_ij_pre,Dn_ij_pre,Dr_ij_projection_pre,Dn_ij_projection_pre]=statistic_P_in_spatical_loc(loc_pre,all_loc);
        [P_post_z,H_post_z,P_post_t,H_post_t,Dr_ij_post,Dn_ij_post,Dr_ij_projection_post,Dn_ij_projection_post]=statistic_P_in_spatical_loc(loc_post,all_loc);
        save([checkpath(fullfile(savepath,'p_spatial\')),[seg_batchi{batchi},'-type',num2str(type)],'-P_spatial.mat'],'loc_pre','loc_post','all_loc','ID_pre','ID_post',...
            'P_pre_z','H_pre_z','P_pre_t','H_pre_t','Dr_ij_pre','Dn_ij_pre','Dr_ij_projection_pre','Dn_ij_projection_pre',...
            'P_post_z','H_post_z','P_post_t','H_post_t','Dr_ij_post','Dn_ij_post','Dr_ij_projection_post','Dn_ij_projection_post', '-v7.3');
    end
end
%% 重新计算p_spatial-types
if true
    for type=1:2
    for batchi=1:3
        for zz=1:length(seg1)
            for hh=1:length(seg2)
                pp=[checkpath(fullfile(savepath,'p_spatial\')),[seg_batchi{batchi},'-type',num2str(type)],'-P_spatial.mat'];
                if exist(pp)
                    load(pp);
                    [P_z_pre,H_z_pre,Dr_ij_projection_pre,Dn_ij_projection_pre]=adjust_P_spatial_loc(loc_pre,Dr_ij_pre,Dn_ij_pre);
                    [P_z_post,H_z_post,Dr_ij_projection_post,Dn_ij_projection_post]=adjust_P_spatial_loc(loc_post,Dr_ij_post,Dn_ij_post);
                    save([checkpath(fullfile(savepath,'p_spatial\')),[seg_batchi{batchi}],'-P_spatial_adjust_projectioned.mat'],'loc_pre','loc_post','all_loc','ID_pre','ID_post',...
                        'P_pre_z','H_pre_z','P_pre_t','H_pre_t','Dr_ij_pre','Dn_ij_pre','Dr_ij_projection_pre','Dn_ij_projection_pre',...
                        'P_post_z','H_post_z','P_post_t','H_post_t','Dr_ij_post','Dn_ij_post','Dr_ij_projection_post','Dn_ij_projection_post', '-v7.3');
                end
            end
        end
    end
    end
end
%% plot pre_post CR spatial_loc
for type=1:2
    for pp=[0.05,0.025]
        for batchi=1:3
            load([checkpath(fullfile(savepath,'p_spatial\')),[seg_batchi{batchi},'-type',num2str(type)],'-P_spatial.mat']);
            for ss=1:2
                switch ss
                    case 1
                        P=P_pre_z;ID=ID_pre;P_r=reshape(P,1,[]);ID_r=reshape(ID,1,[]);
                        P_r(find(ID_r==0))=[];ID_r(find(ID_r==0))=[];
                        l1=['-pre Cond. CR-p'];l2=['-pre Cond. CR'];
                        supervoxel= double(CR_ind_all_pre{batchi,type}(ID_r,1:3));
                    case 2
                        P=P_post_z;ID=ID_post;P_r=reshape(P,1,[]);ID_r=reshape(ID,1,[]);
                        P_r(find(ID_r==0))=[];ID_r(find(ID_r==0))=[];
                        l1=['-post Cond. CR-p'];l2=['-post Cond. CR'];
                        supervoxel= double(CR_ind_all_post{batchi,type}(ID_r,1:3));
                end
                ind=find(P_r<pp);
                temp_plot( supervoxel(ind,:),temp_env.env.supervoxel(:,1:3),temp_env.env.vol,[seg_batchi{batchi},l1,num2str(pp),'-type',num2str(type)],1,checkpath(fullfile(savepath,'p_spatial')),'2021',[])
                
                In=(1-P_r);In=floor(In*100);
                temp_plot( supervoxel,temp_env.env.supervoxel(:,1:3),temp_env.env.vol,[seg_batchi{batchi},l2,'-type',num2str(type)],2,checkpath(fullfile(savepath,'p_spatial')),'2021',In);
            
            end
        end
    end
end


