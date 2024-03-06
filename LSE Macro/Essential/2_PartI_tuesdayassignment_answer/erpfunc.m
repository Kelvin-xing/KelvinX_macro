function erp = erpfunc(coef,par,grid,gh,k)

%Allocating memory
p = zeros(grid.size,1);
er = zeros(grid.size,1);
rf = zeros(grid.size,1);
temp = zeros(gh.size,1);

%Expected rate of return
for i = 1:grid.size,
    p(i,1) = pfunc(grid.d(i,1),coef);
    for j = 1:gh.size,
        dnew = par.mud+par.rhod*grid.d(i,1)+sqrt(2)*par.sigma*gh.e(j,1);
        pnew = pfunc(dnew,coef);
        temp(j,1) = gh.w(j,1)/sqrt(pi)*(dnew+pnew);
    end
    er(i,1) = sum(temp)/p(i,1)-1;
end

%Risk-free rate
for i = 1:grid.size,
    for j = 1:gh.size,
        dnew = par.mud+par.rhod*grid.d(i,1)+sqrt(2)*par.sigma*gh.e(j,1);
        qnew = par.beta*(dnew/grid.d(i,1))^-par.gamma;
        temp(j,1) = gh.w(j,1)/sqrt(pi)*qnew;
    end
    rf(i,1) = 1/(sum(temp))-1;
end


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