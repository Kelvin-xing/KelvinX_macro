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
c^(-nu) = beta*c(+1)^(-nu)*(alpha*z(+1)*k^(alpha-1) + 1 - delta);
c + k   = y + k(-1)*(1-delta);
y       = z*k(-1)^alpha;
z       = 1 - rhoz + rhoz*z(-1) + e;
end;

initval;
k   = ((1-beta*(1-delta))/(beta*alpha))^(1/(alpha-1));
c   = k^alpha - delta*k;
y   = k^alpha;
z   = 1;
end;

shocks;
var e; stderr sigz;
end;

resid;
steady;

estimated_params;
stderr e, 0.01, 0, 0.2;
rhoz, 0.95, 0, 1;
alpha, 0.36, 0, 1;
delta, 0.025, 0, 0.2;
end;

varobs y;
estimation(datafile=y,mode_check) c, k, y;