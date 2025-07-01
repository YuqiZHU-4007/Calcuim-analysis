x = size(env.vol,1);
y = size(env.vol,2);
z=size(env.vol,3);
bg = uint16(100*ones(x,y,z));    %set the background threshould accroding to the background
vol_cut_bg = env.vol - bg;    %substract the background 
figure;imshow(vol_cut_bg(:,:,1),[0 300]);    %test the image quality
for ii = 1 : 3
    aa = medfilt2(vol_cut_bg(:,:,ii),[4 4]);    %wipe off the salt/pepper noise 
    vol_cut_bg_non_salt(:,:,ii) = aa; 
end 
figure;imshow(vol_cut_bg_non_salt(:,:,1),[0 300]);    %test the image quality after wipe off the noise
env.vol(:,:,1:3) = vol_cut_bg_non_salt;   %redefine the env.vol

% vol_signal = uint16(zeros(x,y,z));
env.vol(1:2048,1:685,:) = 0;
env.vol(1:2048,1350:2048,:) = 0;
env.vol(1690:2048,1:2048,:) = 0;
env.vol(160:573,685:1400,:) = 0;
% env.vol(1060:1300,685:1400,:) = 0;
% bb = env.vol.*vol_signal;
% env.vol = bb;
figure;imshow(env.vol(:,:,10),[0 50]);
