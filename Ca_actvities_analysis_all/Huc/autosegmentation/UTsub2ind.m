function k = UTsub2ind(i, j, n)
% UpperTrianlgeIndex(i, j, n)
% n is the number of rows or columns of the square matrix
%

    if i < j
        k = (2*n-i).*(i-1)/2 + (j-i);
    else
        k = (2*n-j).*(j-1)/2 + (i-j);
    end
end

