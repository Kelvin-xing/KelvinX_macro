function [ff] = find_lognormal_difference(sigma,omega,Fomega)
% this function just return the difference of the two sides
% of lognormal distrbution when given sigma, omega, and Fomega

ff=logncdf(omega,-sigma^2/2,sigma)-Fomega;