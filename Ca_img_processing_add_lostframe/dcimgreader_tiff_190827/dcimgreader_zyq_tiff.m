classdef dcimgreader_zyq_tiff < handle %Image4DReader
    % create a dcimgreader that can read from a dcimg file as a 4D (x-y-z-t) image
    % 2016-05-24
    
    properties
        fhandle
        T
        t0
        nz
        nt
        nframes
        height
        width
        T_lessframe
        T_moreframe
        z_lessframe
        z_moreframe
        
        subsamplez
        subsamplet
    end
    
    methods
        function obj = dcimgreader_zyq_tiff(imgpath, T, t0, process ,filetypeF)
            % obj = dcimgreader(imgpath, T, t0 = 1)
            % create a dcimgreader that can read from a dcimg file as a 4D (x-y-z-t) image
            %   t0 = 0 means auto estimate t0
            
            if nargin < 3
                t0 = 0;  
            end
            % inputs
            list=dir([imgpath,'/**/*.tif']);
            obj.fhandle = imfinfo(fullfile([list(1).folder '/' list(1).name]));
            obj.T = T;
            obj.t0 = t0;
            
            obj.T_lessframe=process.T_lessframe;
            obj.T_moreframe=process.T_moreframe;
            obj.z_lessframe=process.z_lessframe;
            obj.z_moreframe=process.z_moreframe;
            % common
            obj.height = obj.fhandle.Height;
            obj.width = obj.fhandle.Width;
            obj.fhandle.imgpath=imgpath;
            obj.fhandle.filetypeF=filetypeF;
            obj.nframes = length(list)-getnonXnum(obj.z_moreframe,0)+getnonXnum(obj.z_lessframe,0);
            

            
            if obj.T == 1  % single plane recording
                if isempty(t0) || t0 == 0
                    obj.t0 = 1;
                end
                obj.nz = 1;
                obj.nt = obj.nframes;
            else  % 3D recording
                if isempty(t0) || t0 == 0  % auto estimate t0
                    obj = estimate_t0(obj);
                end
                % init parameters
                obj.nt = floor(double(obj.nframes - obj.t0 + 1) / obj.T);
                obj.nz = obj.T;
            end
            %
            obj.subsamplez = 1:obj.nz;
            obj.subsamplet = 1:obj.nt;
        end
        
        function obj = estimate_t0(obj)  % only need inputs
            vol = zeros(obj.fhandle.Height, obj.fhandle.Width, obj.T+1, 'uint16');
            for ii = 1:obj.T+1
                frame = tiffread_zyq_0827(obj, ii);
                vol(:,:,ii) = frame;
            end
            totalint = squeeze(sum(sum(vol, 1), 2));%��ά��ȥ��ά��Ϊ1��
            [~, backind] = max(totalint);
            obj.t0 = backind + 1;
            obj.nt = floor(double(obj.nframes - obj.t0 + 1) / obj.T); %
        end
        
        function frame = readframe(obj, ti, zi)%%%�ɸ�
            mf=0;mz=0;lf=0;lz=0;
            indm=find(obj.T_moreframe==ti);indl=find(obj.T_lessframe==ti);
            if isempty(obj.T_lessframe) && isempty(obj.T_moreframe)%%%�治���ڶ�֡��֡���?
                fi = obj.t0 + obj.T * (ti-1) + zi - 1;
                frame = tiffread_zyq_0827(obj, fi);
            else %���ڶ�֡��֡
                if (~isempty(indl) && ~isempty(find(obj.z_lessframe(indl,:)==zi))) %%%ti��zi������,����֡�����?
                    frame = NaN(obj.fhandle.Height, obj.fhandle.Width);%, 'uint16'
                else %����
                    indmz=find(obj.T_moreframe<ti);indlz=find(obj.T_lessframe<ti);
                    if isempty(indmz)
                        mf=0;
                    else
                        mf=getnonXnum(obj.z_moreframe(indmz,:),0);
                        if isempty(indm)
                            mz=0;
                        else
                            mz=getnonXnum(obj.z_moreframe(indm,find(obj.z_moreframe(indm,:)<zi)),0);
                        end
                    end
                    if isempty(indlz)
                        lf=0;
                    else
                        lf=getnonXnum(obj.z_lessframe(indlz,:),0);
                        if isempty(indl)
                            lz=0;
                        else
                            lz=getnonXnum(obj.z_lessframe(indl,find(obj.z_lessframe(indl,:)<zi)),0);
                        end
                    end
                    fi=mf+mz-lf-lz+obj.t0 + obj.T * (ti-1) + zi - 1;
                    fii = obj.t0 + obj.T * (ti-1) + zi - 1;
                    frame = tiffread_zyq_0827(obj, fi);
                end
            end
        end
        
        function vol = readvolume(obj, ti)
%             indl=find(obj.T_lessframe==ti);indm=find(obj.T_moreframe==ti);
%             if length(indl)==0 && length(indm)==0
%                obj.set_subsamplez;
%             elseif length(indl)~=0 && length(indm)==0
%                     list=(1:obj.T);
%                     for i=1:length(indl)
%                         list(find(list-indl(i)==0))=[];
%                     end
%                     obj.set_subsamplet(list);
%             end
            vol = zeros(obj.height, obj.width, numel(obj.subsamplez), 'uint16');
            for ii = 1:numel(obj.subsamplez)%number of elements
                zi = obj.subsamplez(ii);
                frame = obj.readframe(ti, zi);
                vol(:,:,ii) = frame;
            end
        end
        
        function vol = average_volume(obj)
            vol = zeros(obj.height, obj.width, numel(obj.subsamplez));
            for ti = 1:numel(obj.subsamplet)
                a=obj.readvolume(obj, ti);a(:,:,find(isnan(a(1,1,:))))=0;
                vol = vol + a;
            end
            vol = uint16(vol);
        end
        
        function tstack = readslice(obj, zi)
            tstack = zeros(obj.height, obj.width, numel(obj.subsamplet), 'uint16');
            for ii = 1:numel(obj.subsamplet)
                ti = obj.subsamplet(ii);
                frame = obj.readframe(ti, zi);
                tstack(:,:,ii) = frame;
            end
        end
        
        function im4d = read4d(obj)
            
        end
        
        function set_subsamplez(obj, subsamplez)
            if nargin < 2
                subsamplez = 1:obj.nz;
            end
            obj.subsamplez = subsamplez;
        end
        
        function set_subsamplet(obj, subsamplet)
            if nargin < 2
                subsamplet = 1:obj.nt;
            end
            obj.subsamplet = subsamplet;
        end
        
        function close(obj)
            delete(obj);%dcimgclose(obj.fhandle);
        end
    end
    
end








