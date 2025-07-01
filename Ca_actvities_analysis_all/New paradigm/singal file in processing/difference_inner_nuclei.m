%核团内部：左右
[inputname,inputpath]=uigetfile('G:\data\.xls','location');
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
[actname,actpath]=uigetfile('G:\data\.mat','activities_aft_process');
load([behaviorpath,behaviorname]);
load([actpath,actname]);
outputpath=uigetdir('G:\data','outputpath');

[fs_ca,fs_behavior,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);


[~,sheet,~]=xlsfinfo([inputpath,inputname]);sheet
activities_all_region={};
for kk=1:size(sheet,2)
    disp(sheet{kk});
    outpath=checkpath([outputpath '\' sheet{kk}]);
    num=xlsread([inputpath,inputname],sheet{kk});
    num(isnan(num))=0;num(num(:,4)==1,:)=[];
    if isempty(num)
        continue;warning([sheet{kk} 'is empty']);
    end
    %     acti_left=acti(:,find(num(:,3)==1));%左
    %     acti_right=acti(:,find(num(:,3)~=1));%右
    %核团内左右
    checkpath([outpath '\' 'inner nuclei']);
    layer=unique(num(:,1));layer=sort(layer); 
    ind_frames=1:trial.total*frame.per_cycle; %trial.acq(3)*frame.per_cycle+1:trial.test(3)*frame.per_cycle;
    %trial.spon_bef(3)*frame.per_cycle+1:trial.hab(3)*frame.per_cycle; 
    %trial.hab(3)*frame.per_cycle+1:trial.acq(3)*frame.per_cycle;
    y_left=[];
    for jj=1:2
       switch jj
            case 1
                ind_left=num(find(num(:,3)==1),2);%左
                layer_left=num(find(num(:,3)==1),1);
                lab='left';
            case 2
                ind_left=num(find(num(:,3)~=1),2);%右
                layer_left=num(find(num(:,3)~=1),1);
                lab='right';
        end
        h_left=figure;title([sheet{kk} '-' lab]);
        set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.88]) ;
        %ind_right=num(find(num(:,1)==layer(ii) & num(:,3)~=1),:);%右  
        for ii=1:length(layer_left)
            y_left{jj}(:,ii)=activities_preCS_dfdf_aftcorrect{layer_left(ii),1}(ind_frames,ind_left(ii));
        end
        str=strcat({char(string(layer_left))},{char(char('-')*ones(size(ind_left)))},{char(string(ind_left))});
        [~,~,sig]=sepplot_20181222(ind_frames,y_left{jj},str);hold on;%title([sheet{kk} '-' num2str(layer(ii)) 'layer-left'])
        %US line
        line([trial.hab(3)*frame.per_cycle+frame.us_start:frame.per_cycle:trial.acq(3)*frame.per_cycle;trial.hab(3)*frame.per_cycle+frame.us_start:frame.per_cycle:trial.acq(3)*frame.per_cycle]',...
            [min(sig(:)) max(sig(:))],'color','r','linestyle','--');hold on;
        %CS bar
        patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
            [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
            [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
            [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
            repmat([min(sig(:)) min(sig(:)) max(sig(:)) max(sig(:))],trial.test(3)-trial.hab(2)+1,1)',...
            'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.15);hold on
        xlim([ind_frames(1) ind_frames(end)]);
        saveas(h_left,[outpath '\inner nuclei\' sheet{kk} '-' lab '.tif']);%close(h_left);
    end
    
    h_mean=figure;
    %US line
    line([trial.hab(3)*frame.per_cycle+frame.us_start:frame.per_cycle:trial.acq(3)*frame.per_cycle;trial.hab(3)*frame.per_cycle+frame.us_start:frame.per_cycle:trial.acq(3)*frame.per_cycle]',...
        [min(y_left{1}(:)) max(y_left{1}(:))],'color','r','linestyle','--');hold on;
    %CS bar
    patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
        repmat([min(y_left{1}(:)) min(y_left{1}(:)) max(y_left{1}(:)) max(y_left{1}(:))],trial.test(3)-trial.hab(2)+1,1)',...
        'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.15);hold on
    p1=plot(ind_frames,mean(y_left{1},2),'r');hold on;unnamed(:,2)=mean(y_left{1},2);
    p2=plot(ind_frames,mean(y_left{2},2),'b');
    legend([p1,p2],{'left','right'});
    xlim([ind_frames(1) ind_frames(end)]);
        
    %核团间均值和onset
    acti=[];a_event={};a_event_rep=[];
    activities_all_region{kk,1}=sheet{kk};
    num_smallsize=num(:,find(num(:,4)~=1));
    if ~isempty( num_smallsize)
        for ii=1:size(num_smallsize,1)
            ind_frames=trial.hab(3)*frame.per_cycle+1:trial,acq(3)*frame.per_cycle;
            acti(:,ii)=activities_preCS_dfdf_aftcorrect{ num_smallsize(ii,1),1}(ind_frames, num_smallsize(ii,2));
            activities_all_region{kk,2}=mean(acti,2);
            a_event{ii,1}=activities_event_preCS{ num_smallsize(ii,1),1}{ num_smallsize(ii,2),1}';
            a_event_onset1=[];
            a_event_onset2=[];
            onset=[];
            onset(:,1)=ceil(a_event{ii,1}.ind/frame.per_cycle);
            onset(:,2)=mod(a_event{ii,1}.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
            a_event_onset1=[a_event_onset1;onset(:,1)];
            a_event_onset2=[a_event_onset2;onset(:,2)];
        end
        a_event_onset=[a_event_onset1 a_event_onset2];summ=[];
        for jj=trial.hab(3)+1:trial.acq(3)
            ind=find(a_event_onset(:,1)==jj & a_event_onset(:,2)<frame.cs_end & a_event_onset(:,2)>=frame.cs_start);
            if ~isempty(ind) &&  mean(a_event_onset(ind,2))<frame.us_start
                summ.m(jj-trial.hab(3))=mean(a_event_onset(ind,2));
                summ.sem(jj-trial.hab(3))=std((a_event_onset(ind,2)-frame.cs_start));%/sqrt(length(ind));
                summ.ind(jj-trial.hab(3))=true;
                summ.mx(jj-trial.hab(3))=jj;
            else
                summ.m(jj-trial.hab(3))=0; summ.ind(jj-trial.hab(3))=0;summ.mx(jj-trial.hab(3))=jj;summ.sem(jj-trial.hab(3))=0;
            end
        end
        activities_all_region{kk,3}=[summ.ind;summ.m;summ.sem]';
    else
        warning(['there are no singel neuron inner ' sheet{kk}]);
    end
end

%不同核团均值，acq
figure,
for kk=1:size(sheet,2)
    plot( activities_all_region{kk,2});hold on;
end

for ii=1:trial.acq(1)
    figure,
    for kk=1:size(sheet,2)
        plot(activities_all_region{kk,2}(1+(ii-1)*frame.per_cycle:ii*frame.per_cycle));hold on
        scatter(activities_all_region{kk,3}(ii,2),activities_all_region{kk,2}(activities_all_region{kk,3}(ii,2)));hold on;
        line([activities_all_region{kk,3}(ii,2)-activities_all_region{kk,3}(ii,3) activities_all_region{kk,3}(ii,2)+activities_all_region{kk,3}(ii,3)],...
            [activities_all_region{kk,2}(activities_all_region{kk,3}(ii,2)) activities_all_region{kk,2}(activities_all_region{kk,3}(ii,2))]);hold on
    end
end
