function xyrot = imrotate_xy(im, ang, xy)
% xyrot = imrotate_xy(im, ang, xy)
% crop mode is needed

    ori = (size(im) + 1) / 2;
    ori = ori(end:-1:1);
    rotmat = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
    oriall = repmat(ori, size(xy,1), 1);
    xyrot = (xy - oriall) * rotmat + oriall;
end
