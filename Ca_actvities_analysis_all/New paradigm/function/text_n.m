function text_n(n,B,B_control_CS,B_control_US,figure_able,figure_able_contolCS,figure_able_contolUS)
a=0.1;pos=max(max(B(:))-a*min(B(:)),0);
pos_CS=max(max(B_control_CS(:))-a*min(B_control_CS(:)),0);
pos_US=max(max(B_control_US(:))-a*min(B_control_US(:)),0);
posX=size(B,1)-0.5;
if figure_able
    if ~isempty(B)
        B=cur_nan_col(B')';
        if  sum(n==[8:9])==2 || sum(n==3)==1
            text(posX,pos,['n=' num2str(size(B,2)/2)],'fontsize',20,'interpreter','none');
        else
            text(posX,pos,['n=' num2str(size(B,2)/2)],'fontsize',20,'interpreter','none');
        end
    end
end
if  figure_able_contolCS
    if ~isempty(B_control_CS)
        B_control_CS=cur_nan_col(B_control_CS')';
        if  sum(n==[8:9])==2 || sum(n==3)==1
            text(posX,pos-3*a*pos,['nCS=' num2str(size(B_control_CS,2)/2)],'fontsize',20,'interpreter','none');
        else
            text(posX,pos-3*a*pos,['nCS=' num2str(size(B_control_CS,2)/2)],'fontsize',20,'interpreter','none');
        end
    end
end
if figure_able_contolUS
    if ~isempty(B_control_US)
        B_control_US=cur_nan_col(B_control_US')';
        if  sum(n==[8:9])==2 || sum(n==3)==1
            text(posX,pos-6*a*pos,['nUS=' num2str(size(B_control_US,2)/2)],'fontsize',20,'interpreter','none');
        else
            text(posX,pos-6*a*pos,['nUS=' num2str(size(B_control_US,2)/2)],'fontsize',20,'interpreter','none');
        end
    end
end