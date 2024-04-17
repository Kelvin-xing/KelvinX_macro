% A code for implementing the HP filter

function [tt,ct] = hp_filter(y,lambda)

[T,n] = size(y);

LAM = zeros(T,T);

LAM(1,1) = 1+lambda;
LAM(1,2) = -2*lambda;
LAM(1,3) = lambda;

LAM(2,1) = -2*lambda;
LAM(2,2) = 1+ 5*lambda;
LAM(2,3) = -4*lambda;
LAM(2,4) = lambda;

for j = 3:T-2
    LAM(j,j) = 1 + 6*lambda;
    LAM(j,j+1) = -4*lambda;
    LAM(j,j+2) = lambda;
    LAM(j,j-1) = -4*lambda;
    LAM(j,j-2) = lambda;
end

LAM(T,T) = 1 + lambda;
LAM(T,T-1) = -2*lambda;
LAM(T,T-2) = lambda;
LAM(T-1,T) = -2*lambda;
LAM(T-1,T-1) = 1 + 5*lambda;
LAM(T-1,T-2) = -4*lambda;
LAM(T-1,T-3) = lambda;

%tt = ((LAM)^(-1))*y;
tt = LAM\y;
%tt = inv(LAM)*y;
ct = y - tt;

