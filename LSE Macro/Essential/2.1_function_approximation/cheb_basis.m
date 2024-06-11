function xout = cheb_basis(x,N)

%  
%  constructs first N Chebyshev basis functions (including zeroth order)
%

dim = size(x,1);
xout = ones(dim,N);
xout(:,2) = x;

for n = 3:N+1
    xout(:,n) = 2*x.*xout(:,n-1)-xout(:,n-2);
end
