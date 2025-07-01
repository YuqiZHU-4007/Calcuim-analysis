function number=getnonXnum(A,X)
if isempty(A)
    number=0;
else
number=numel(A(find(A~=X)));
end
end