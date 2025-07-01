function [A,ind]=adjust_location(AA,res,adjust_matrix)
ind=find(AA(:,2)<adjust_matrix(3) & (floor(AA(:,3)/res(3)+1)>=adjust_matrix(1) &floor(AA(:,3)/res(3)+1)<=adjust_matrix(2)));
A=AA;
A(ind,2)=AA(ind,2)+adjust_matrix(4);
end

% seg{1}={'20220803fish2','20220814fish1','20210402fish1'};%L
% seg{2}={'20210702fish1'};%C
% seg{3}={'20220824fish1','20220829fish1','20210529fish1'};%NL
% seg{4}={'20210717fish1','20210822fish1'};%FL
% adjust_matrix{1}=[[14,31,570,40];[12,22,525,40];[14,22,525,55]];
% adjust_matrix{2}=[18,22,500,30];
% adjust_matrix{3}=[[15,21,575,40];[15,19,540,40];[13,18,540,35]];
% adjust_matrix{4}=[[15,21,530,45];[16,26,540,50]];
% 
% columns = {'x', 'y', 'z', 't','label'};
% res=[0.66,0.66,10];
% slicea=reg_mask; sliceb = rescalegd(im2double(slicea(:,:,1:end-1)), [1/10000 1/10000]);
% outpath=checkpath('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\');
for batchi=1:length(seg)
    for jj=1:length(seg{batchi})
        nn=seg{batchi}{jj};
        AA=readmatrix(['H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP\',nn,'vol_env_spatialloc_warped_SyN.csv']);
        outpathh=[outpath,[nn,'vol_env_spatialloc_warped_SyN.csv']]
        [AA_aft,ind]=adjust_location(AA,res,adjust_matrix{batchi}(jj,:));
        M=table(AA_aft(:,1),AA_aft(:,2),AA_aft(:,3),AA_aft(:,4),AA_aft(:,5),'VariableNames', columns);
        writetable(M,outpathh);
        %check plot
        A=[floor(AA(:,1)/res(2)),floor(AA(:,2)/res(1)),floor(AA(:,3)/res(3)+1)];%raw
        B=[floor(AA_aft(:,1)/res(2)),floor(AA_aft(:,2)/res(1)),floor(AA_aft(:,3)/res(3)+1)];%aft
        showspv = zeros(size(slicea,1), size(slicea,2), 3, size(sliceb,3));
        showspv_aft = zeros(size(slicea,1), size(slicea,2), 3, size(sliceb,3));
        for ii=1:size(sliceb,3)
            pti=find(A(:,3)==ii);
            slice = insertShape(sliceb(:,:,ii), 'circle', [abs(A(pti, 1)), abs(A(pti, 2)), 5*ones(length(pti),1)]);
            showspv(:,:,:,ii) = slice;%showspv(:,:,2,ii) = slice;showspv(:,:,3,ii) = slice;
             pti=find(B(:,3)==ii);
            slice = insertShape(sliceb(:,:,ii), 'circle', [abs(B(pti, 1)), abs(B(pti, 2)), 5*ones(length(pti),1)]);
            showspv_aft(:,:,:,ii) = slice;%showspv(:,:,2,ii) = slice;showspv(:,:,3,ii) = slice;
        end
        show_spv_GUI(showspv);title(nn)
        show_spv_GUI(showspv_aft);title(nn)
    end
end


