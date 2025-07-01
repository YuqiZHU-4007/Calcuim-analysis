function idx = remove_overlap(pos, rad, scores)
%function idx = remove_overlap(points, fim, radmap)
%
% TODO: not proper input

%{
[yy, xx] = find(points);
np = length(xx);

rad = radmap(:,:,1);
rad = rad(points);
scores = fim(points);
[~, order] = sort(scores, 'descend');

pos = [xx, yy];
rad = rad;
scores;
%}

np = size(pos, 1);

k = min([10 np]);

[~, order] = sortrows([-scores, abs(rad - mean(rad))]);

[nearest, dis] = knnsearch(pos, pos, 'k', k);
%nearest = nearest(:,2:end);
%dis = dis(:,2:end);

idx = zeros(np, 1);
for ii = 1:np
    pi = order(ii);
    for jj = 2:k
        ni = nearest(pi,jj);
        d = dis(pi,jj);
        if rad(ni) + rad(pi) >= d  && any(idx == ni)
            break;
        end
    end
    if jj == k
        idx(ii) = pi;
    end
end
idx = sort(nonzeros(idx));



        
    