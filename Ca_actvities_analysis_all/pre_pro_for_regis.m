%pre process for registration
cut_move=false;
T=1:24;
img_pth={}
img_pth{1}='/mnt/X/data/fear_conditioning/huc/20190905/fish1/spon_1/rec168087208/raw tiff/';
img_pth{3}='/mnt/X/data/fear_conditioning/huc/20190905/fish1/spon_2/rec172554109/raw tiff/';
img_pth{2}='//mnt/X/data/fear_conditioning/huc/20190905/fish1/rec169059079/raw tiff/';

img_out_pth='/mnt/X/data/fear_conditioning/huc/20190905/fish1/for regis/';
beh_path_1='/mnt/X/data/fear_conditioning/huc/20190905/fish1/spon_1/shake.mat';
beh_path_trainig='/mnt/X/data/fear_conditioning/huc/20190905/fish1/20190905_fish1_huc.mat';
beh_path_2='/mnt/X/data/fear_conditioning/huc/20190905/fish1/spon_2/shake.mat';
para_path='/mnt/X/data/fear_conditioning/huc/20190905/fish1/para.mat';
dsample_meth='box';
%load('E:\A_Data_lightsheet\Data_huc\20190905\fish1\training\para.mat');
% stimCS=zeros(1,frameb.per_cycle*trial.total(1));
% for ii=1:size(stimCS,2)/frameb.per_cycle
%     stimCS((ii-1)*frameb.per_cycle+frameb.cs_start:(ii)*frameb.per_cycle)=0;%23
% end
%% load
ref={};
figure,
for ii=1:3
    refi=[];
    switch ii
        case 1
            ref_num=1:10;
            load(beh_path_1);
        case 2
            ref_num=1:10;
            load(beh_path_trainig);
        case 3
            ref_num=1:10;
            load(beh_path_2);
    end
    a=load(para_path);
    y=false(size(y_3sd));y(y_3sd~=0)=true;y(y_3sd==0)=false;
    if ii==2
        for jj=1:3
            switch jj
                case 1
                    frame_ind=(1:a.trial.hab(3)*a.frameb.per_cycle);
                case 2
                    frame_ind=(a.trial.hab(3)*a.frameb.per_cycle+1:a.trial.acq(3)*a.frameb.per_cycle);
                case 3
                    frame_ind=(a.trial.acq(3)*a.frameb.per_cycle+1:a.trial.test(3)*a.frameb.per_cycle);
            end
            b=find(~y(frame_ind));b=b(1:ceil(a.fs.behavior):length(b));refi=[refi,frame_ind(b(ref_num))];
        end
    else
        b=find(~y);b=b(1:ceil(a.fs.behavior):length(b));refi=b(ref_num);
    end
    subplot(3,1,ii),plot(y_3sd);hold on;
    line([refi;refi]',[-5,5],'color','k','linewidth',3);xlim([1,length(y_3sd)]);
    ref{ii}=ceil(refi/a.fs.behavior/a.fs.ca);
end
%% downsample

%% save as nrrd
I={};I_d={};
for ii=1:3
    ind=ref{ii};
    if ~isfolder(fullfile([img_out_pth,'session_',num2str(ii)]))
        mkdir(fullfile([img_out_pth,'session_',num2str(ii)]));
    end
    for jj=1:length(ind)
        for t=T
            I{ii,jj}(:,:,t)=imread(fullfile([img_pth{ii},'z',num2str(t,'%02d'),'/',num2str(ind(jj),'%04d'),'.tif']));
            I_d{ii,jj}(:,:,t)=imresize(I{ii,jj}(:,:,t), [1348/4 NaN],dsample_meth);
%             figure,imshow(I{ii,jj}(:,:,t),[0 500]);
%             figure,imshow(I_d{ii,jj}(:,:,t),[0 500]);
        end
        %     tif=Tiff([img_out_pth,'session_',num2str(ii),'\',num2str(ind(jj),'%04d'),'.Tiff'],'w');
        %     setTag(tif,tagstruct)
        %     write(tif,I{ii,jj});
        saveastiff(I_d{ii,jj}, fullfile([img_out_pth,'session_',num2str(ii),'/',num2str(ind(jj),'%04d'),'.tif']))
    end
end
save([img_out_pth 'I.mat'],'I','ref','I_d','-v7.3');

% I=imread([img_pth,'z',num2str(t,'%02d'),'\',num2str(ind(jj),'%04d'),'.tif']);
% figure,
% imshow(I,[0 500]);
% met={'nearest','bilinear','bicubic','lanczos2','lanczos3','box','pyramid'};
% for ii=1:size(met,2)
%     if strcmp(met{ii},'pyramid')
%         I_d=impyramid(I, 'reduce');
%         I_d=impyramid(I_d, 'reduce');
%     else
%         I_d=imresize(I, [1348/4 NaN],met{ii});
%     end
%     figure,imshow(I_d,[0 500]);title(met{ii});
% end