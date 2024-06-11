function xout = pol_basis(x,N)

%  
%  constructs first N regular polynomial basis functions
%

dim = size(x,1);
xout = ones(dim,N);

for n = 2:N+1
    xout(:,n) = x.^(n-1);
end
