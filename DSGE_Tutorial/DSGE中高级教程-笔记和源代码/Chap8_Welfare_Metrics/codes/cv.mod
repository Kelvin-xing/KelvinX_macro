/*This file try to implement the ideas in the Chap8 Welfare Metric
 *All variables are in levels; we totally have 12 endogenous variables 
 *The utility form taken from the example in section 2.3.2, additively separable
*/
var w   //welfare metric
w_c    //consumption part
w_l   // labor
k      //capital stock
a     //technology shock
c     //consumption
n    //labor
wage   //real wage 
Rk   //capital rate
y    //output
i    //investment
r   //real rate
;

varexo e;

parameters alpha beta delta rho sigma psi phi sigmae;
parameters ws wcs wls ks as cs ns wages Rks ys is rs;

alpha = 1/3; % capital share
beta = 0.995; %discout factor
delta = 0.02; %capital depreciation rate
rho =0.95; % persistence of technology shock
sigma = 1.05; %CRRA parameter
phi = 0.4; % inverse of labor Frisch elasitcity
%sigmae = 0.01; % lower standard deviation of technology shock;

load parameterfile_cv;
set_param_value('sigmae',sigmae);

as = 1;
ns = 1/3;
ks = ns*(alpha/(1/beta - 1 +delta))^(1/(1-alpha));
is = delta*ks;
Rks = as*alpha*(ks)^(alpha-1)*ns^(1-alpha);
wages =as*(1-alpha)*(ks)^alpha*ns^(-alpha);
ys  =as*(ks)^alpha*ns^(1-alpha);
cs = ys - is;
psi= wages*cs^(-sigma)/ns^phi ;
rs = Rks - delta;
wcs =  (cs^(1-sigma) - 1)/(1-sigma)/(1-beta);
wls = - psi*ns^(1+phi)/(1+phi)/(1-beta);
ws = wcs  + wls;
model;
%(1) marginal product of capital stock = capital rate
Rk = a*alpha*(k(-1))^(alpha-1)*n^(1-alpha);

%(2)marginal product of labor = real wage
wage = a*(1-alpha)*(k(-1))^alpha*n^(-alpha);

%(3) technology shock
a = (1-rho)+rho*a(-1) + e;

%(4)capital accumulation
i = k - (1-delta)*k(-1);

%(5)labor supply equation
psi*n^phi = wage*c^(-sigma);

%(6) output
y = a*(k(-1))^alpha*n^(1-alpha);

%(7)resource constraint
y = c+ i;

%(8) Euler equation
c^(-sigma)  = c(+1)^(-sigma)*beta*(1+r);

%(9) real rate
1+r  = Rk(+1) + 1- delta;

%(10) welfare metric
w = w_c + w_l;

%(11) consumption part
w_c = (c^(1-sigma) - 1)/(1-sigma) + beta*w_c(+1);

%(12) labor part
w_l = - psi*n^(1+phi)/(1+phi) + beta*w_l(+1);
end;

initval;
w = ws;
w_c = wcs;
w_l = wls;
k = ks ;
a =as;
c = cs;
n = ns;
wage = wages;
y = ys;
i = is;
Rk = Rks;
r = rs;
end;

shocks;
var e = sigmae^2;
end;

steady;
resid(1);

stoch_simul(order =2, irf=0);