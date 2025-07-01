function [ROI_traces_all,ROI_ind] = extra_ROI_of_selectregion(extra_signal,num_type)
% automatically chose ROIs
% 201411_YangChen
%20141201_cyin
%20190102_zyq : Write index of component for definded regions
extra_signal=0;
num_type='layered';
%%import vol
[fn_tif,fn_path] = uigetfile('E:\A_Data_lightsheet\Data_vmat/*.tiff','Vol'); % using tif stack file
%%import mat
[env_name,env_path] = uigetfile([fn_path '/*.mat'],'env','MultiSelect', 'off'); % using tif stack file
[act_name,act_path] = uigetfile([fn_path '/*.mat'],'act','MultiSelect', 'off'); % using tif stack file

%load act
if extra_signal
    [act_name,act_path] = uigetfile([fn_path '/*.mat'],'act','MultiSelect', 'off'); % using tif stack file
    selectfile=fullfile(act_path,act_name);
    for ii=1:size(selectfile,2)
        load(selectfile{ii});
    end
end
%outpath
outpath=checkpath([fn_path,'subregoin_supervol']);
load([env_path env_name]);
%read stack
stack = my_tiffread([fn_path,fn_tif]);%imread(fn_tif); %my_tiffread.m same as imread.m for tiff stack frame by frame
%%  parameters
para.region_num=input('input the number of regions you selected to extra signal\n');
supervoxel=env.supervoxel;
region_name=cell(para.region_num,1);
outpath_region=cell(para.region_num,1);
edge_x=cell(para.region_num,1);edge_y=cell(para.region_num,1);
mask_region=cell(para.region_num,1);
for zz=1:para.region_num
    region_name{zz}=input('input the name of regions you selected to extra signal\n','s');
    para.T{zz}=input(['input the layer of ', region_name{zz},' regions, eg. 2:5\n']);
    outpath_region{zz}=checkpath([outpath '\' region_name{zz}]);
    %manually define ROI area
    for ii= para.T{zz}
        slice =stack(:,:,ii);
        slice = normim(slice);
        gs = fspecial('gaussian', opt.maxrad, 1);
        slice = imfilter(slice, gs);
        
        h1=figure('name',['layer' num2str(ii) '-' region_name{zz}]); imshow(slice);
        l=true;
        for ll=1:2
        if l
            [BW1,c1,r1]= roipoly(slice);
            l=false;
        else 
            slice(BW1)=min(slice(:));
            [BW2,c2,r2]= roipoly(slice);
        end
        end
        BW= BW1 | BW2;
        c={c1 c2};r={r1 r2};
        if ~isempty(BW)
            mask_region{zz}(:,:,ii) = BW; % manually define ROI area
            edge_x{zz}{ii,1}=c;edge_y{zz}{ii,1}=r;
        else
            mask_region{zz}(:,:,ii)=zeros(size(slice));
            edge_x{zz}{ii,1}=0;edge_y{zz}{ii,1}=0;
        end
        h2=figure; imshow(mask_region{zz}(:,:,ii));
        close(h1);close(h2);
    end
end
save([fn_path '/' 'ROI_infor.mat'], 'mask_region','edge_x','edge_y','region_name', '-v7.3');

%%select supervolex
ROI_trace_all=cell(para.region_num,1);
for zz=1:para.region_num
    ROI_ind=cell(length(para.T{zz}),1); ROI_traces=cell(length(para.T{zz}),1);
    ROI_ind_rep=[];
    for ii= para.T{zz}
        spvind=find(supervoxel(:,3) == ii)';
        spv=supervoxel(spvind,1:2);
        ind1=find(mask_region{zz}(:,:,ii));
        ind2 = sub2ind([env.height env.width], spv(:,2), spv(:,1));
        [ind,i1.i2]=intersect(ind1,ind2);
        [ind_inregion_y,ind_inregion_x]=ind2sub([env.height env.width],ind);
        if ~isempty(ind_inregion_x)
            for jj=1:length(ind_inregion_x)
                ROI_ind{ii,1}(jj)=find(supervoxel(:,1)==ind_inregion_x(jj) & supervoxel(:,2)==ind_inregion_y(jj) & supervoxel(:,3)==ii);
                switch num_type
                    case 'raw'
                        ROI_ind_rep=[ROI_ind_rep;ii ROI_ind{ii,1}(jj)'];
                    case 'layered'
                        a=ROI_ind{ii,1}(jj)'-(length(find(supervoxel(:,3)<ii)));
                        ROI_ind_rep=[ROI_ind_rep;ii a];
                end
                if extra_signal
                    ROI_traces{ii,1}(:,jj)=activities(ROI_ind{ii,1}(jj),:);
                end
            end
        end
        f=figure('Name','autoROI','color', 'white', 'position', [400   200   size(slice)]);% 1 plot
        slice =stack(:,:,ii);
        slice = normim(slice);
        slice(ind)=true;
        imshow(slice);hold on;
        text(ind_inregion_x, ind_inregion_y, string(num2str((1:length(ind_inregion_x))')),'fontsize',11,'HorizontalAlignment','center');hold on;%ROI_ind{ii,1}(:)
        %plot(edge_x{zz}{ii,1},edge_y{zz}{ii,1},'r','linewidth',1.5);hold on;
        rad=(supervoxel(ROI_ind{ii,1},4)+supervoxel(ROI_ind{ii,1},5))/2;
        ellipse(rad,rad,zeros(size(rad)),ind_inregion_x,ind_inregion_y,1.5);
        saveas(f,[outpath_region{zz},'\',num2str(ii),'.tif']);
        saveas(f,[outpath_region{zz},'\',num2str(ii),'.fig']);
        close(f);
        %save([fn_path '/' 'ROI_traces_' num2str(zz) '.mat'], 'ROI_traces', 'ROI_ind','region_name', '-v7.3');
    end
    if extra_signal
        ROI_trace_all{zz,1}=ROI_traces; ROI_trace_all{zz,2}=ROI_ind; ROI_trace_all{zz,3}=region_name{zz};
    end
    %
    xlswrite([outpath '/location'],ROI_ind_rep,region_name{zz});
end
if extra_signal
    save([fn_path '/' 'ROI_trace_all.mat'], 'ROI_trace_all', '-v7.3');
else
    repmat_act_to_cell_format([act_path act_name],[env_path env_name]);
end

