function [C,sumd,D] = FindCentroid_Direct(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
D = zeros(numU,1);
for i = 1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
        sumd(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,sum,D1] = kmeans(M_s,1,'distance','sqeuclidean');
        C(i,:) = C1;
        D(i) = mean(D1);
        sumd(i)=sum;
    end
end
end