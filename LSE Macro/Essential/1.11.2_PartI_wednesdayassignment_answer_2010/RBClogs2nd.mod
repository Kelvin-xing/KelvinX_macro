var c, k, z;

varexo e;

parameters alpha, beta, delta, nu, rho, sigma;

alpha = 0.33;
beta  = 0.96;
delta = 0.04;
rho   = 0.98;

load nuparam;
set_param_value('nu',nu)
load sigparam;
set_param_value('sigma',sigma)


model;
exp(-nu*c) = beta*exp(-nu*c(+1))*(alpha*z(+1)*exp((alpha-1)*k) + 1 - delta);
exp(c) + exp(k) = z*exp(alpha*k(-1)) + (1-delta)*exp(k(-1));
z = (1-rho) + rho*z(-1) + e;

end;

initval;
k = log(((1/beta-1+delta)/alpha)^(1/(alpha-1)));
c = log(exp(alpha*k) - delta*exp(k));
z = 1;
 
end;

steady;

check;

shocks;
var e; stderr sigma;
end;

stoch_simul(order=2,nocorr,nomoments,nograph,IRF=40);