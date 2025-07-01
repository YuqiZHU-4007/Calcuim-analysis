function frame = readframe_zyq(imio, ti, zi)
% frame = readframe(imio, ti, zi)
%
% do not support vector index so far
%

switch imio.filetype
    case 'dcimg',
%         obj.T_lessframe=opt.process.T_lessframe;
%         obj.T_moreframe=opt.process.T_moreframe;
%         obj.z_lessframe=opt.process.z_lessframe;
%         obj.z_moreframe=opt.process.z_moreframe;
        indm=find(imio.T_moreframe==ti);indl=find(imio.T_lessframe==ti);
        if isempty(imio.T_lessframe,2) && isempty(imio.T_moreframe,2)%%%存不存在多帧少帧情况
            fi = imio.t0 + imio.T * (ti-1) + zi - 1;
            frame = dcimgmatlab(fi-1, imio.basepath);
        else %存在多帧少帧
            if (length(indl)~=0 && length(find(imio.z_lessframe(indl,:)==zi))~=0) %%%ti的zi不存在,在少帧情况里
                frame = NaN(imio.height, imio.width);%, 'uint16'
            else %存在
                indmz=find(imio.T_moreframe<ti);indlz=find(imio.T_lessframe<ti);
                if isempty(indmz)
                    mf=0;
                else
                    mf=getnonXnum(imio.z_moreframe(indmz,:),0);
                    if isempty(indm)
                        mz=0;
                    else
                        mz=getnonXnum(imio.z_moreframe(indm,find(imio.z_moreframe(indm,:)<zi)),0);
                    end
                end
                if isempty(indlz)
                    lf=0;
                else
                    lf=getnonXnum(imio.z_lessframe(indlz,:),0);
                    if isempty(indl)
                        lz=0;
                    else
                        lz=getnonXnum(imio.z_lessframe(indl,find(imio.z_lessframe(indl,:)<zi)),0);
                    end
                end
                fi=mf+mz-lf-lz+imio.t0 + imio.T * (ti-1) + zi - 1;
                frame = dcimgmatlab(fi-1, imio.basepath);
            end
        end
        %fi = imio.t0 + imio.T * (ti-1) + zi - 1;
        
        frame = frame';
        clear functions;
        
    case 'mptiff',
        fi = imio.t0 + imio.T * (ti-1) + zi - 1;
        frame = imread(imio.basepath, fi);
        
    case 'tif',
        fi = imio.t0 + imio.T * (ti-1) + zi - 1;
        frame = imread(fullfile(imio.basepath, imio.readpath{fi}));
        
    case 'vol',
        frame = imread(fullfile(imio.basepath, imio.readpath{ti}), zi);
        
    case 'zt',
        frame = imread(fullfile(imio.basepath, imio.readpath{ti, zi}));
end

end



