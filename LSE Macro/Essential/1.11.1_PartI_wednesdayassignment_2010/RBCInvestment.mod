var c, k, inv, y, z;

varexo e;

parameters alpha, beta, delta, nu, rho, sigma;

alpha = 0.33;
beta  = 0.96;
delta = 0.04;
nu    = 2;
rho   = 0.98;
sigma = 0.007;


model;
exp(-nu*c)/exp(z) = beta*exp(-nu*c(+1))*(alpha*exp((alpha-1)*k) + (1 - delta)/exp(z(+1)));
exp(c) + exp(inv) = exp(alpha*k(-1));
exp(k)            = (1-delta)*exp(k(-1)) + exp(z)*exp(inv);
exp(y)            = exp(alpha*k(-1));
z                 = rho*z(-1) + e;

end;

initval;
k   = log(((1/beta-1+delta)/alpha)^(1/(alpha-1)));
c   = log(exp(alpha*k) - delta*exp(k));
inv = log(delta*exp(k));
y   = log(exp(alpha*k));
z   = 0;
 
end;

steady;

check;

shocks;
var e; stderr sigma;
end;

stoch_simul(order=1,nocorr,nomoments,IRF=40);