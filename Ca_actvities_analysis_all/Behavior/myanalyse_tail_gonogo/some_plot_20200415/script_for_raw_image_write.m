clc;clear all
addpath(genpath('/mnt/Z/Behavior/myanalyse_tail_gonogo/some_plot_20200415/plot_CS_US_co'))
[fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=0;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num)
%plot_behavior_onset(a.delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
%spon
T_spon={};T_US={};l_spon={};l_US={};T_non_spon=[];
for ii=1:trial.acq_block_num+2
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

[n,p]=uigetfile('/mnt/w/data/fear conditioning_ZYQ/*.seq');
% 31,32,33,34,35,36

T=[32,];
for t=unique(T)
    b1=checkpath(fullfile(p,[n,'raw_Image_Trail',num2str(t,'%02d')]));
    for ii=(t-1)*frameb.per_cycle+1 : t*frameb.per_cycle
        [y, ~] = readtailseq([p,'/',n],T_non_spon(ii));
        if ii>= (t-1)*1800+721 && ii <= (t-1)*1800+1009
            y=insertShape(y,'filledcircle',[160 225 30],'color','r','Opacity',1);
        end
        %figure,imshow(y)
        imwrite(y,[b1,'/',num2str(ii,'%05d'),'.tif']);
    end
end

%
% p=uigetdir('/mnt/X/data/fear_conditioning/huc/');
% % n=fileparts(p);
% %31,32,33,34,35,36
% 
% T=[1,2,3,4,5,6];
% for t=unique(T)
%     b1=checkpath(fullfile([n,'/raw_Image_Trail',num2str(t,'%05d')]));
%     for ii=(t-1)*frameb.per_cycle+1 : t*frameb.per_cycle
%         [y, ~] = imread([p,'/',num2str(T_non_spon(ii),'%05d'),'.tif']);
%         if ii>= (t-1)*frameb.per_cycle+frameb.cs_start && ii <= (t-1)*frameb.per_cycle+frameb.cs_end
%             y=insertShape(y,'filledcircle',[160 225 30],'color','r','Opacity',1);
%         end
%         %figure,imshow(y)
%         imwrite(y,[b1,'/',num2str(ii,'%05d'),'.tif']);
%     end
% end