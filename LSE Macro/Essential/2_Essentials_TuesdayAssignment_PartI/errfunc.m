function [ssr error] = errfunc(coef,par,grid,gh)
% note that the second argument is ignored by optimizer

%QUESTION 1.1. COMPLETE THIS FUNCTION.

%Allocating memory
lhs = zeros(grid.size,1);
rhs = zeros(grid.size,1);
temp = zeros(gh.size,1);

%Double for loop
for i = 1:grid.size,
    %CALCULATE HERE THE LHS OF THE EULER EQUATION. HINT: YOU CAN USE PFUNC.
    %lhs(i,1) = ...;
    for j = 1:gh.size,
        %CALCULATE HERE THE RHS OF THE EULER EQUATION FOR EACH GAUSSIAN
        %HERMITE NODE.
        %temp(j,1) = ...;
    end
    %CALCULATE HERE THE RHS OF THE EULER EQUATION BY COMBINING THE ELEMENTS
    %OF TEMP. (NOTE THAT THE ELEMENTS OF TEMP WILL BE OVERWRITTEN BY THE
    %CALCULATIONS FOR THE NEXT GRID POINT.)
    %rhs(i,1) = ...;
end

%Sum of squared Euler equation errors
%CALCULATE HERE THE SUM OF SQUARED EULER EQUATION ERRORS.
%ssr = ...;

error = lhs-rhs;

end