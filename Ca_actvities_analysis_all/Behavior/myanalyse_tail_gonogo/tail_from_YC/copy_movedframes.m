

imdir = 'Z:\data\20151110\20151110_tail';
cpdir = checkpath('Z:\data\20151110\tail_moved');
thres = 4.2e-8;


moving = (stdwintail > thres);
movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));
movind = find(movlabel);

figure; plot(stdwintail);
hold on; plot(moving * thres, 'r');

for fi = movind(:)'
    fi
    filename = [num2str(fi-1) '.bmp'];
    if exist([cpdir '/' filename], 'file') > 0
        continue;
    end
    copyfile([imdir '/' filename], [cpdir '/' filename]);
end
    