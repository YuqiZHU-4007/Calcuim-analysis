function [i, j] = UTind2sub(k, N)
% [i, j] = UTind2sub(k, N)
% N is the length of the 1-d representation of the matrix
%

    n = (sqrt(8*N+1) + 1) / 2;

    b = - 2*n - 1;
    ac = 2*n + 2*k;

    i = ceil((-b - sqrt(b^2 - 4*ac)) / 2 - 1);

    j = k - (2*n-i).*(i-1)/2 + i;

end

