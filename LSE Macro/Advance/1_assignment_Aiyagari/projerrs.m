function ssr = projerrs(coef,par,grid,gh,siz)

% left hand side of Euler equation (depending on order of polynomial)
c = polyn(grid.a,coef,par.order);

%Right hand side of Euler equation
rhs = zeros(siz.a,siz.e);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below fill in the necessary steps to obtain the rhs of the Euler equation
% for each Gauss-Hermite node and for all points in the state space
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for i = 1:siz.e
    
    knew = ;
    anew = ;
    % tomorrows consumption (depending on order of polynomial)
    cnew = polyn();
    rhs(:,i) = ;
end
rhss = sum(rhs,2);

%Norm
ssr = norm(c.^(-par.nu)-rhss);

end