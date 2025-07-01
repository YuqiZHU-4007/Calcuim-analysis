function US_mapping_singletrial_fraction(Path)
set(0,'defaultfigurecolor','w')
addpath(genpath('F:\DUlab\FC_analyse\FishExplorer'));
if isempty(Path)
    p=uigetdir('H:\1.Test US\2.Tail free¡ª¡ªData from 117\')
else
    p=Path
end
if (exist([fullfile(p,'singletrial'),'\CR_ind_summary_singletrial.mat'],'file')==2)
    load(fullfile(p,'para.mat'));
    load(fullfile(p,'env.mat'));
    load(fullfile(p,'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
    nn=[Path(end-14:end-7),Path(end-5:end-1)];
    p=checkpath(fullfile(p,'singletrial'));
    %warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
     warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
    supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
%    clr=jet(31);a=str2double(brain_region_id(:,1));b=isnan(a);a(b)=31;
%     figure,scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,[0.5 0.5 0.5],'filled');hold on;
%     for ii=1:length(index_in_region_in_clust);
%         ind=index_in_region_in_clust{ii};
%         scatter3(supervolxeli(ind,1),supervolxeli(ind,2),supervolxeli(ind,3),10,clr(a(ind),:),'filled');hold on;
%     end
    
    load(fullfile(p,'CR_ind_summary_singletrial.mat'), 'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS_inthissession');
    ind_all_emergeCSUS=cell(length(labels_typei),trial.total,trial.total,2,1);
    
    Fraction_in_region_type=cell(length(labels_typei),2,2);Index_in_region=cell(length(labels_typei),trial.total,2,2);
    Loc_in_region_cell=cell(length(labels_typei),trial.total,2,2);Num_in_region=cell(length(labels_typei),2,2,3);
    Brain_region_id=cell(length(labels_typei),trial.total,2,2);Label_region=cell(length(labels_typei),2,2);
    
    Fraction_in_region_type_emerged=cell(length(labels_typei),2,2,1);Index_in_region_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);
    Loc_in_region_cell_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);Num_in_region_emerged=cell(length(labels_typei),2,2,3,1);
    Brain_region_id_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);Label_region_emerged=cell(length(labels_typei),2,2,1);
    
    Fraction_in_region_type_emerged_inthissession=cell(length(labels_typei),2,2,3);Index_in_region_emerged_inthissession=cell(length(labels_typei),trial.total,2,3);
    Loc_in_region_cell_emerged_inthissession=cell(length(labels_typei),trial.total,2,2,3);Num_in_region_emerged_inthissession=cell(length(labels_typei),2,2,3,3);
    Brain_region_id_emerged_inthissession=cell(length(labels_typei),trial.total,2,2,3);Label_region_emerged_inthissession=cell(length(labels_typei),2,2,3);
    
    sessionx = {'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'};
    for iscutmov=1
        %% ÄÔÇøÍ³¼Æ
        for jj=1:length(labels_typei)
            for ii=1:1:trial.total
                for zz=1:trial.total
                    a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                    b=ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov};
                    ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}=setdiff(a,b);
                end
            end
        end
        %fraction
        fraction_type=1;
        for jj=1:length(labels_typei)
            for ii=1:trial.total
                ind_type=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                Brain_region_id{jj,ii,iscutmov,fraction_type}=brain_region_id(ind_type,:);
                for regioni=1:length(Label)
                    if length(index_in_region_in_clust{regioni,1})<=10
                        ind_inregion=[];
                    else
                        ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                    end
                    Index_in_region{jj,ii,iscutmov,fraction_type}{regioni,:}=ind_inregion;
                    Fraction_in_region_type{jj,iscutmov,fraction_type}(regioni,ii)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                    Loc_in_region_cell{jj,ii,iscutmov,fraction_type}{regioni,:}=supervolxeli(ind_inregion,:);
                    Num_in_region{jj,iscutmov,fraction_type,1}(regioni,ii)= length(ind_inregion)./length(ind_type);
                    Num_in_region{jj,iscutmov,fraction_type,2}(regioni,ii)= length(ind_inregion);
                    Num_in_region{jj,iscutmov,fraction_type,3}(regioni,ii)= length(ind_type);
                    Label_region{jj,iscutmov,fraction_type}(regioni,:)=string(Label{regioni});
                end
            end
        end
        for jj=1:length(labels_typei);%[1:18,31:34]
            for typei=1:size(ind_all_emergeCSUS_inthissession,4)
                for ii=1:trial.total
                    ind_type=ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,typei};
                    Brain_region_id_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}=brain_region_id(ind_type,:);
                    for regioni=1:length(Label)
                        if length(index_in_region_in_clust{regioni,1})<=10
                            ind_inregion=[];
                        else
                            ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                        end
                        Index_in_region_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}{regioni,:}=ind_inregion;
                        Fraction_in_region_type_emerged_inthissession{jj,iscutmov,fraction_type,typei}(regioni,ii)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                        Loc_in_region_cell_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}{regioni,:}=supervolxeli(ind_inregion,:);
                        Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,1,typei}(regioni,ii)= length(ind_inregion)./length(ind_type);
                        Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,2,typei}(regioni,ii)= length(ind_inregion);
                        Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,3,typei}(regioni,ii)= length(ind_type);
                        Label_region_emerged_inthissession{jj,iscutmov,fraction_type,typei}(regioni,:)=string(Label{regioni});
                    end
                end
            end
        end
        for zz=[]
            for jj=1:length(labels_typei);%[1:18,31:34]
                for  typei=1:size(ind_all_emergeCSUS,5)
                    for ii=1:trial.total
                        for zz=1:trial.total
                            ind_type=ind_all_emergeCSUS{jj,ii,zz,iscutmov,typei};
                            Brain_region_id_emerged{jj,ii,zz,iscutmov,fraction_type,typei}=brain_region_id(ind_type,:);
                            for regioni=1:length(Label)
                                if length(index_in_region_in_clust{regioni,1})<=10
                                    ind_inregion=[];
                                else
                                    ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                                end
                                Index_in_region_emerged{jj,ii,zz,iscutmov,fraction_type,typei}{regioni,:}=ind_inregion;
                                Fraction_in_region_type_emerged{jj,iscutmov,fraction_type,typei}(regioni,ii,zz)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                                Loc_in_region_cell_emerged{jj,ii,zz,iscutmov,fraction_type,typei}{regioni,:}=supervolxeli(ind_inregion,:);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,1,typei}(regioni,ii,zz)= length(ind_inregion)./length(ind_type);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,2,typei}(regioni,ii,zz)= length(ind_inregion);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,3,typei}(regioni,ii,zz)= length(ind_type);
                                Label_region_emerged{jj,iscutmov,fraction_type,typei}(regioni,:)=string(Label{regioni});
                            end
                        end
                    end
                end
            end
        end
        save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'ind_all_emergeCSUS','Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region',...
            'Fraction_in_region_type_emerged_inthissession','Index_in_region_emerged_inthissession','Loc_in_region_cell_emerged_inthissession','Num_in_region_emerged_inthissession','Brain_region_id_emerged_inthissession','Label_region_emerged_inthissession',...
            'Fraction_in_region_type_emerged','Index_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged',...
            '-append','-v7.3');
        
    end
end
end
