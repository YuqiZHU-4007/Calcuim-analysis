imio = dcimgreader_zyq(opt.imgpath, opt.T, [],opt.process);  %, opt.t0); % auto estimate t
savepath = '';
poolsize = ;
%delete(gcp('nocreate'));
poolobj = parpool(poolsize);

tic;
spmd(poolsize)
    for zi = labindex:poolsize:env.depth
        [labindex zi]
        zpath = checkpath([savepath '/raw tiff/z' num2str(zi, '%02d')]);
        for ti = 1:env.nvol
            frame = readframe(imio, ti, zi);
            imwrite(uint16(frame), [zpath '/' num2str(ti, '%04d') '.tif']);
        end
    end
end
imio.close();
toc;