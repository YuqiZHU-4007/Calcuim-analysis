clear all;
path= uigetdir('/scratch/20191214/')%'E:\A_Data_lightsheet\Data_huc\20190905\fish2\brain area\';
[actname,actpath]=uigetfile(path,'env');load([actpath actname]);
[ filelist] = scanDir_types(path,'.txt');
reg_mask=zeros(env.height,env.width,env.depth);reg_name={};reg_loc=struct;
for jj=1:length(filelist)
    [~,~,x,y,z] = textread(filelist{jj},'%s%s%s%s%s');%,'headerlines',1
    [~,reg_name{jj}]=fileparts(filelist{jj});
    xx=nan(1,length(x));yy=nan(1,length(x));zz=nan(1,length(x));
    %l=cell2table(l,'VariableNames',{'x','y','z' });
    for ii=1:length(x)
        a=(x{ii});
        if ~isempty(a)
            if  ~strcmp(x{ii},'X') && ~strcmp(x{ii},'') && str2num(a)>0 && str2num(y{ii})>0
                xx(ii)=str2num(a);
                yy(ii)=str2num(y{ii});
                zz(ii)=str2num(z{ii});
                reg_mask(xx(ii),yy(ii),zz(ii))=jj;
            else
                disp([reg_name{jj} ' ' num2str(ii) '_' x{ii} '_' y{ii} '_' z{ii}]); continue; 
            end
        else
            continue;
        end
    end
    xx(isnan(xx))=[];    yy(isnan(yy))=[];zz(isnan(zz))=[];
    if isfield(reg_loc,strrep(reg_name{jj},'-','_'))
        reg_loc=setfield(reg_loc,strrep(reg_name{jj},'-','_'),[]);
    end
    reg_loc=setfield(reg_loc,strrep(reg_name{jj},'-','_'),[xx;yy;zz]);
end
for jj=1:length(filelist)
    reg_name{jj}=strrep(reg_name{jj},'_',' ');
end
%reg_mask(env.width+1:end,:,:)=[];reg_mask(:,env.height+1:end,:)=[];
save([path, 'region_mask.mat'],'reg_mask','reg_loc','reg_name');
show_spv_GUI(reg_mask);