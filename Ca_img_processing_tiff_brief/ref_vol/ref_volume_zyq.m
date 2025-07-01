function vol = ref_volume_zyq(imio, nsamples, zi,process,savepath)
% vol = ref_volume(imio, nsamples, zi)
% 
% generate reference volume using average of volumes
% who did not move much according to zi slice.
%

if nargin < 3
    zi = ceil(imio.T / 2);
end
if nargin < 2
    nsamples = 128;
end

% % tstack = imio.readslice(zi);
% if length(opt.T2)~=0
%     list=(1:nsamples);
%     for i=1:length(opt.T2)
%     list(find(list-opt.T2(i)==0))=[];
%     end
%     imio.set_subsamplet(list);
% else
%     imio.set_subsamplet(1:nsamples);
% end
%             obj.T_lessframe=process.T_lessframe;
%             obj.T_moreframe=process.T_moreframe;
%             obj.z_lessframe=process.z_lessframe;
%             obj.z_moreframe=process.z_moreframe;


list=newsubsamplet(nsamples,zi,process);
 imio.set_subsamplet(list);
 tstack = imio.readslice(zi);

% imio.read4d();
% tic; tstack = dr.readslice(17); toc;  % caution

% save cache for checking reasons
seqwrite(tstack, [savepath, '/z' num2str(zi, '%.2d') '_for_vol/'])

difff = zeros(size(tstack,3), 1);
for ii = 2:size(tstack,3)
%     difff(ii) = sum(sum(tstack(:,:,ii) - tstack(:,:,ii-1)));
    sig1 = double(tstack(:,:,ii));
    sig0 = double(tstack(:,:,ii-1));
    difff(ii) = corr(sig1(:), sig0(:));%每张与前一张的corr(一个值)
end
    
h=figure('visible','off'); plot(difff(2:end));
saveas(h,[savepath '/corr1.tif']);
% if max(difff(2:end))<=0.8
%     error('bad reference in nsample and zi');
%     displayresult('bad reference in nsample and zi');
% end

sv = sort(difff);

isstation = difff > sv(round(nsamples/4));
% 75% top correlation to previous frame

temp = bwlabel(isstation);
temp2 = hist(temp, unique(temp));
temp3 = find(temp2 <= 3) - 1;
for ii = temp3(:)'
    isstation(temp==ii) = 0;%把与其他不在一个连通区域并且个数较少的层去掉
end

stationind = find(isstation);%Find the nonzero elements

meanslice = mean(tstack(:,:,stationind), 3);%%得到的corr较大的层数取mean

diffm = zeros(size(tstack,3), 1);
for ii = 1:size(tstack,3)
%     difff(ii) = sum(sum(tstack(:,:,ii) - tstack(:,:,ii-1)));
    sig1 = double(tstack(:,:,ii));
    sig0 = meanslice;
    diffm(ii) = corr(sig1(:), sig0(:));
end
h=figure('visible','off'); plot(diffm(stationind))
saveas(h,[savepath '/corr2.tif']);

[sv, order] = sort(diffm);
selti = intersect(order(ceil(nsamples/2):end), stationind);%Create two vectors that have some values in common
% % top correlation to meanslice

if isempty(selti)
    warning('bad reference');%%%判断选出来的几张与这几张和mean的50%的corr有没有重叠
    displayresult('bad reference in average stack');
    selti = stationind;
end

% tic;
vol = double(imio.readvolume(selti(1)));
for ii = selti(2:end)'
    vol = vol + double(imio.readvolume(ii));
end
% toc;
vol = vol / length(selti);%%计算所有层的average mask
vol = uint16(vol);
% % % ImageJ;
% % % MIJ.createImage(vol);

%tiffwrite(uint16(vol), [opt.savepath '/vol.tiff']);









