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
c^(-nu) = beta*c(+1)^(-nu)*(alpha*z(+1)*k^(alpha-1) + 1 - delta);
c + k = z*k(-1)^alpha + (1-delta)*k(-1);
z = (1-rho) + rho*z(-1) + e;

end;

initval;
k = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
c = k^alpha - delta*k;
z = 1;
 
end;

steady;

check;

shocks;
var e; stderr sigma;
end;

stoch_simul(order=1,nocorr,nomoments,nograph,IRF=40);