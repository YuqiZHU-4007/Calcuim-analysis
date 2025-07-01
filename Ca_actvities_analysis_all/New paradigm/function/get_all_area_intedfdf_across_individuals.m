name=fieldnames(CS_acq_block)
Acq_intedfdf_dim_flash_m=[];Acq_intedfdf_dim_flash_sd=[];Acq_intedfdf_dim_flash_n=[];Acq_intedfdf_dim_flash_area_raw=[];
Acq_intedfdf_dim_flash_m_hab5=[];Acq_intedfdf_dim_flash_sd_hab5=[];Acq_intedfdf_dim_flash_n_hab5=[];Acq_intedfdf_dim_flash_area_raw_hab5=[];
hab_tst_intedfdf_dim_flash_m=[];hab_tst_intedfdf_dim_flash_sd=[];hab_tst_intedfdf_dim_flash_n=[];hab_tst_intedfdf_dim_flash_area_raw=[];
m1=[];sd1=[];m=[];sd=[];n=[];
for ii=1:length(name)
    ind=getfield(CS_acq_ind,name{ii})';ind=unique(ind);
    
    areaa=getfield(CS_acq_block,name{ii})';
    m=[];sd=[];n=[];
    for jj=ind
        if jj>length(actpathes_dim)
            win=1:10;
        else
            win=1:5;
        end
        ind_filed=find(getfield(CS_acq_ind,name{ii})==jj);
        areaa_filed=areaa(win,ind_filed);
        m(win,jj)=mean(areaa_filed,2);
        sd(win,jj)=std(areaa_filed,[],2);
        n(jj,1)=size(areaa_filed,2);
        area_raw{jj,1}=areaa_filed;
    end
    Acq_intedfdf_dim_flash_m=setfield(Acq_intedfdf_dim_flash_m,name{ii},m);
    Acq_intedfdf_dim_flash_sd=setfield(Acq_intedfdf_dim_flash_sd,name{ii},sd);
    Acq_intedfdf_dim_flash_n=setfield(Acq_intedfdf_dim_flash_n,name{ii},n);
    Acq_intedfdf_dim_flash_area_raw=setfield(Acq_intedfdf_dim_flash_area_raw,name{ii},area_raw);
    
    areaa=getfield(CS_hab_tst5,name{ii})';
    m1=[];sd1=[];n1=[];
    for jj=ind
        win=1:2;
        ind_filed=find(getfield(CS_acq_ind,name{ii})==jj);
        areaa_filed=areaa(win,ind_filed);
        m1(win,jj)=mean(areaa_filed,2);
        sd1(win,jj)=std(areaa_filed,[],2);
        n1(win,jj)=size(areaa_filed,2);
        area_raw1{jj,1}=areaa_filed;
    end
    Acq_intedfdf_dim_flash_m_hab5=setfield(Acq_intedfdf_dim_flash_m_hab5,name{ii},[m1(1,:);m]);
    Acq_intedfdf_dim_flash_sd_hab5=setfield(Acq_intedfdf_dim_flash_sd_hab5,name{ii},[sd1(1,:);sd]);
    Acq_intedfdf_dim_flash_n_hab5=setfield(Acq_intedfdf_dim_flash_n_hab5,name{ii},n);
    Acq_intedfdf_dim_flash_area_raw_hab5=setfield(Acq_intedfdf_dim_flash_area_raw_hab5,name{ii},area_raw1);
    
    
    areaa=getfield(CS_hab_tst,name{ii})';
    m=[];sd=[];n=[];
    for jj=ind
        win=1:2;
        ind_filed=find(getfield(CS_acq_ind,name{ii})==jj);
        areaa_filed=areaa(win,ind_filed);
        m(win,jj)=mean(areaa_filed,2);
        sd(win,jj)=std(areaa_filed,[],2);
        n(win,jj)=size(areaa_filed,2);
        area_raw{jj,1}=areaa_filed;
    end
    hab_tst_intedfdf_dim_flash_m=setfield(hab_tst_intedfdf_dim_flash_m,name{ii},m);
    hab_tst_intedfdf_dim_flash_sd=setfield(hab_tst_intedfdf_dim_flash_sd,name{ii},sd);
    hab_tst_intedfdf_dim_flash_n=setfield(hab_tst_intedfdf_dim_flash_n,name{ii},n);
    hab_tst_intedfdf_dim_flash_area_raw=setfield(hab_tst_intedfdf_dim_flash_area_raw,name{ii},area_raw);
end


path='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328';
save([path '\Intedfdf_all_learner.mat'],'Acq_intedfdf_dim_flash_m','Acq_intedfdf_dim_flash_sd','Acq_intedfdf_dim_flash_n','Acq_intedfdf_dim_flash_area_raw',...
    'Acq_intedfdf_dim_flash_m_hab5', 'Acq_intedfdf_dim_flash_sd_hab5', 'Acq_intedfdf_dim_flash_n_hab5', 'Acq_intedfdf_dim_flash_area_raw_hab5',...
    'hab_tst_intedfdf_dim_flash_m', 'hab_tst_intedfdf_dim_flash_sd', 'hab_tst_intedfdf_dim_flash_n', 'hab_tst_intedfdf_dim_flash_area_raw');