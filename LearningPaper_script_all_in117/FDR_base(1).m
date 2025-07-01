function [H,q]=FDR_base(pval,method)
if strcmp(method,'BH2')
    [sorted_p,sort_index]=sort(pval);
    for i=1:length(pval)
        q(i)=pval(i)*length(pval)/sort_index(i);
    end
    H=(q<=0.05);
end

if strcmp(method,'Bonferroni')
    for i=1:length(pval);
        q(i)=pval(i);
    end
    H=(q<=(0.05/length(pval)));
end

if strcmp(method,'BH')
    [sorted_p,sort_index] = sort(pval);
    m = length(pval);
    for i = 1:m
        q(i) = sorted_p(i) * m / i;
    end
    q = q(sort_index);
    H = (q <= 0.05);
end
end
