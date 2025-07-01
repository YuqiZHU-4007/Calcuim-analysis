classdef dcimgreader2 < handle %Image4DReader
%%use reshape 帧同步信号作为新的索引 
%未改    
    properties
        fhandle
        T
        t0
        nz
        nt
        nframes
        height
        width
        
        subsamplez
        subsamplet
    end
    
    methods
        function obj = dcimgreader2(imgpath, T, t0)
            % obj = dcimgreader(imgpath, T, t0 = 1)
            % create a dcimgreader that can read from a dcimg file as a 4D (x-y-z-t) image
            %   t0 = 0 means auto estimate t0
            
            if nargin < 3
                t0 = 0;  
            end
            if ~libisloaded('tmcamcon')
                loadlibrary('tmcamcon.dll', 'tmcamcon.h');
            end
            
            % inputs
            obj.fhandle = dcimgopen(imgpath);
            obj.T = T;
            obj.t0 = t0;
            % common
            obj.height = obj.fhandle.height;
            obj.width = obj.fhandle.width;
            obj.nframes = obj.fhandle.nframes;
            
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
            vol = zeros(obj.fhandle.height, obj.fhandle.width, obj.T+1, 'uint16');
            for ii = 1:obj.T+1
                frame = dcimgread(obj.fhandle, ii);
                vol(:,:,ii) = frame;
            end
            totalint = squeeze(sum(sum(vol, 1), 2));%降维，去掉维度为1的
            [~, backind] = max(totalint);
            obj.t0 = backind + 1;
            obj.nt = floor(double(obj.nframes - obj.t0 + 1) / obj.T); %
        end
        
        function frame = readframe(obj, ti, zi)
            fi = obj.t0 + obj.T * (ti-1) + zi - 1;
            frame = dcimgread(obj.fhandle, fi);
        end
        
        function vol = readvolume(obj, ti)
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
                vol = vol + obj.readvolume(obj, ti);
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
            dcimgclose(obj.fhandle);
        end
    end
    
end








