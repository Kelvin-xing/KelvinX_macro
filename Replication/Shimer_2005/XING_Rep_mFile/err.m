%%%% This file calculates the error in the fitting procedure in
%%%% Approximation method

%%%% coef is the coefficient vector to optimize over
%%%% F_new is the maximized new value function
%%%% p is the productivity grid

function err = err(coef, F_new, p)
    err = 0;
    [n,~] = size(p);
    for i = 1:n
        approxi = coef(1) + p(i) * (coef(2) + p(i) * (coef(3) + p(i) * (coef(4) + p(i) * coef(5))));
        err = err + (approxi-F_new(i))^2;
    end
end