classdef Image4DReader

    properties
        fhandle
        T
        t0
        nz
        nt
        nframes
        %nBins
        
        subsamplez
        subsamplet
    end

    methods
        function obj = Image4DReader(imgpath, T, t0)
            % switch path
            
            
            
            
            obj.fhandle = dcimgread(imgpath);
            obj.T = T;
            obj.t0 = t0;
            obj.nz = T;
            obj.nframes = obj.fhandle.nframes;
            obj.nt = floor(double(nframes - t0 + 1) / T);
        end
        
        function frame = readframe(ti, zi)
            fi = imio.t0 + imio.T * (ti-1) + zi - 1;
            frame = dcimgreadframe(fhandle, fi);
        end
        
        function vol = readvolume(ti, subsamplez)
            for ii = 1:numel(subsamplez)
                zi = subsamplez(ii);
                frame = readframe(imio, ti, zi);
                vol(:,:,ii) = frame;
            end
        end
        
        function tstack = readslice(zi, subsamplet)
            
        end
        
    end
end



