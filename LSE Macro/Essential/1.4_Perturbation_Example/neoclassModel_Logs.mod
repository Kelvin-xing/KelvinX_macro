// neoclassical growth model solution and simulation

var c, k, y, z;
varexo e;
parameters alpha, beta, delta, nu, rhoz, sigz, kss, css;

load params;

set_param_value('alpha'     ,par.alpha);   // returns to scale parameter
set_param_value('beta'      ,par.beta);    // discount factor
set_param_value('delta'     ,par.delta);   // depreciation rate
set_param_value('nu'        ,par.nu);      // relative risk aversion coefficient
set_param_value('rhoz'      ,par.rhoz);    // autocorrelation of productivity shock
set_param_value('sigz'      ,par.sigz);    // standard deviation of productivity shock

set_param_value('kss'       ,par.k);       // steady state capital
set_param_value('css'       ,par.c);       // steady state consumption

model;
exp(c)^(-nu) = beta*exp(c(+1))^(-nu)*(alpha*exp(z(+1))*exp(k)^(alpha-1) + 1 - delta);
exp(c) + exp(k)   = exp(y) + exp(k(-1))*(1-delta);
exp(y)       = exp(z)*exp(k(-1))^alpha;
exp(z)       = 1 - rhoz + rhoz*exp(z(-1)) + e;
end;

initval;
k   = log(kss);
c   = log(css);
y   = log(exp(kss)^alpha);
z   = 0;
end;

shocks;
var e; stderr sigz;
end;

resid;
steady;

stoch_simul(order=1,nomoments, irf=20, periods = 5000);