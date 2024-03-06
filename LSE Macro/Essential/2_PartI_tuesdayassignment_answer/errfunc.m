function ssr = errfunc(coef,par,grid,gh)

%QUESTION 1.1. COMPLETE THIS FUNCTION.

%Allocating memory
lhs = zeros(grid.size,1);
rhs = zeros(grid.size,1);
temp = zeros(gh.size,1);


%Double for loop
for i = 1:grid.size,
    %CALCULATE HERE THE LHS OF THE EULER EQUATION. HINT: YOU CAN USE PFUNC.
    %lhs(i,1) = ...;
    lhs(i,1) = pfunc(grid.d(i,1),coef);
    for j = 1:gh.size,
        %CALCULATE HERE THE RHS OF THE EULER EQUATION FOR EACH GAUSSIAN
        %HERMITE NODE.
        %temp(j,1) = ...;
        dnew = par.mud+par.rhod*grid.d(i,1)+sqrt(2)*par.sigma*gh.e(j,1);
        pnew = pfunc(dnew,coef);
        qnew = par.beta*(dnew/grid.d(i,1))^-par.gamma;
        temp(j,1) = gh.w(j,1)/sqrt(pi)*qnew*(dnew+pnew);
    end
    %CALCULATE HERE THE RHS OF THE EULER EQUATION BY COMBINING THE ELEMENTS
    %OF TEMP. (NOTE THAT THE ELEMENTS OF TEMP WILL BE OVERWRITTEN BY THE
    %CALCULATIONS FOR THE NEXT GRID POINT.)
    %rhs(i,1) = ...;
    rhs(i,1) = sum(temp);
end

%Sum of squared Euler equation errors
%CALCULATE HERE THE SUM OF SQUARED EULER EQUATION ERRORS.
%ssr = ...;
ssr = norm(lhs-rhs);

end