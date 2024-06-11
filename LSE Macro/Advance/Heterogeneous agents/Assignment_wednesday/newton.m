function [apg,apb] = newton(r,rp,z,zp,ap_g1,ap_b1,ap_g,ap_b,newt_g,newt_b)

global a_grid tau mu T alpha delta beta sigma phi

% This program uses somewhat of a Newton-Raphson method to find the
% solution to the Euler equation today, given some policy function for
% tomorrow.

w = (1-alpha)*exp(z)*((r+delta)/alpha )^(alpha/(alpha-1));

wp = (1-alpha)*exp(zp)*((rp+delta)/alpha )^(alpha/(alpha-1));

metric = 1;

while metric>1e-1

    ap_g1n = newt_g(a_grid,ap_g1,interp1(a_grid,ap_g,ap_g1,'linear','extrap'),max(interp1(a_grid,ap_b,ap_g1,'linear','extrap'),phi),r,w,rp,wp);

    ap_b1n = newt_b(a_grid,ap_b1,interp1(a_grid,ap_g,ap_b1,'linear','extrap'),max(interp1(a_grid,ap_b,ap_b1,'linear','extrap'),phi),r,w,rp,wp);

    metric = max(max(abs([ap_g1-ap_g1n ap_b1-ap_b1n])));

    ap_g1 = ap_g1n; ap_b1 = ap_b1n;

end

apg = ap_g1; apb = ap_b1;

end