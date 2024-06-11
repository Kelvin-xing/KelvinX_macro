function ssr = projerrs_full(coef,par,grid,gh,siz)

% left hand side of Euler equation (depending on order of polynomial)
c = polyn(grid.a,coef,par.order);

%Right hand side of Euler equation
rhs = zeros(siz.a,siz.e);
for i = 1:siz.e
    
    knew = grid.a - c;
    anew = (1 + par.r - par.delta)*knew + par.wss*(1+sqrt(2)*par.sigshock*gh.e(i));
    % tomorrows consumption (depending on order of polynomial)
    cnew = polyn(anew,coef,par.order);
    rhs(:,i) = (gh.w(i)/sqrt(pi))*(-par.zeta2 + par.zeta1*exp(-par.zeta0*knew) ...
        +par.beta*cnew.^(-par.nu)*(par.r + 1 - par.delta));
end
rhss = sum(rhs,2);

%Norm
ssr = norm(c.^(-par.nu)-rhss);

end