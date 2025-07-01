function A=get_integrate_X1_X2(A1,A2)
if ~isempty(A1)
 A=[A1(1,:);A2;A1(2,:)];
else
    A=A2;
end