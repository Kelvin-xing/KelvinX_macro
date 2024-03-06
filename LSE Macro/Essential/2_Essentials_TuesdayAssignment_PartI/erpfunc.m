function erp = erpfunc(coef,par,grid,gh,k)

%Allocating memory
p = zeros(grid.size,1);
er = zeros(grid.size,1);
rf = zeros(grid.size,1);
temp = zeros(gh.size,1);

%Expected rate of return
for i = 1:grid.size,
    p(i,1) = ...
    for j = 1:gh.size,
        ...
    end
    er(i,1) = sum(temp)/p(i,1)-1;
end

%Risk-free rate


if k == 1
    erp = er-rf;
end
if k == 2
    erp = er;
end
if k == 3
    erp = rf;
end


end