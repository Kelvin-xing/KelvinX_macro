// simple model with endogenous labor

var c, k, h, z;
varexo e;

parameters beta, rho, alpha, nu, delta, kappa, init_c, init_k, init_h;

alpha   = 0.33;
rho     = 0.95;
beta    = 0.99;
kappa   = 0.5;
delta   = 0.025;

load parameterfile
set_param_value('nu',nu);
set_param_value('init_c',init_c);
set_param_value('init_k',init_k);
set_param_value('init_h',init_h);

%nu = 0.1;
%nu = 4;
%init_k = 79.6248;
%init_c =  6.47881;
%init_h =  2.80879;

model;
c^(-nu)=beta*c(+1)^(-nu)*(exp(z(+1))*alpha*k^(alpha-1)*h(+1)^(1-alpha) + 1-delta );
c+k=exp(z)*k(-1)^alpha*h^(1-alpha)+(1-delta)*k(-1);
c^(-nu)*(1-alpha)*exp(z)*k(-1)^alpha*h^(-alpha)=h^kappa;
z = rho*z(-1)+e;
end;

initval;
c = init_c;
k = init_k;
h = init_h;
z = 0;
end;

shocks;
var e; stderr 0.007;
end;

steady;

//stoch_simul(order=1,solve_algo=3,nocorr,nomoments,IRF=0);
