function invec = vecin(vec, inset)
% invec = vecin(vec, inset)
%
% return logical vector indicate whether element of vec is in inset
%

invec = arrayfun(@(x) any(x==inset), vec);



