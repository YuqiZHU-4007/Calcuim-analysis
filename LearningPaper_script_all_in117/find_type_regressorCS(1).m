function ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei)
%thrfun=@(x)(0.5);
thrfun1=@(x)(nanmean(x)+2*nanstd(x));%CS
thrfun2=@(x)(nanmean(x)+2*nanstd(x));%US
thrfun3=@(x)(nanmean(x)+3*nanstd(x));%
 
stimcorr=CS_related_stimcorr{batchi,fishi};
thr=thrfun1(stimcorr);
CS_related_ind=find(stimcorr>thr);
 
stimcorr=US_related_stimcorr{batchi,fishi};
thr=thrfun2(stimcorr);
US_related_ind=find(stimcorr>thr);
 
stimcorr=CSUS_related_stimcorr{batchi,fishi};
B=sort( stimcorr,'descend');thr=0.7;B(floor(length(stimcorr)*0.005));
%thr=thrfun2(stimcorr);
CSUS_type=find(stimcorr>thr);
ind_thr=find(sum(Thr_related{batchi,fishi,1}(2:4,:),1)>=2 & sum(Thr_related{batchi,fishi,2}(5:7,:),1)>=2);
 
switch typei
    case 1
        ind_type=CS_related_ind;
    case 2
        ind_type=US_related_ind;
    case 3
        ind_type=CSUS_type;
    case 4
        ind_type=intersect(CS_related_ind,ind_thr);
        %ind_type=intersect(CS_related_ind,US_related_ind);
end
end

