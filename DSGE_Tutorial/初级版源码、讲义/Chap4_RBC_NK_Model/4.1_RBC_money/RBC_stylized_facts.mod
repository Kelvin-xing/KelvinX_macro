%This file is written by Xiangyang Li@Aug.,11,2015
%Simulation to find out success and failures of RBC model.
var a n c k i r w R y yn;
varexo ea;
parameters theta alpha delta beta sigma rho sda;
parameters ns cs ks as is rs ws Rs;

alpha = 0.33;
beta  = 0.99;
delta = 0.025;
sigma = 1;
rho = 0.974; %calibrated
sda = 0.009;  %calibrated

as = 1;
ns = 1/3;
rs = 1/beta;
kn = (alpha/(1/beta - 1 +delta))^(1/(1-alpha));
ks = ns*kn;
Rs= 1/beta - 1 + delta;
ws = (1-alpha)*kn^alpha;
ys =as*kn^alpha*ns;
is=  delta*ks;
cs = ys - is;
theta = (1-ns)/cs*(1-alpha)*as*kn^alpha;

model;
% (1) Euler equation, capital
exp(c)^(-sigma)=beta*exp(c(+1))^(-sigma)*(R(+1)+(1-delta));

% (2) Euler equation, bonds
exp(c)^(-sigma)=beta*exp(r)*exp(c(+1))^(-sigma);

% (3) Labor supply
theta/(1-exp(n))=exp(c)^(-sigma)*exp(w);

% (4) Production func
exp(y)=exp(a)*exp(k(-1))^(alpha)*exp(n)^(1-alpha);

% (5) Capital demand
R=alpha*exp(a)*exp(k(-1))^(alpha-1)*exp(n)^(1-alpha);

% (6) Labor demand
exp(w)=(1-alpha)*exp(a)*exp(k(-1))^(alpha)*exp(n)^(-alpha);

% (7) Resource constraint
exp(y)=exp(c)+exp(i);

% (8) Capital accumulation
exp(k)=exp(i)+(1-delta)*exp(k(-1));

% (9) Productivity shock (TFP)
a=rho*a(-1)+ea;

%(10) labor average productivity
exp(yn) = exp(y)/exp(n);

end;
initval;
k=log(ks);
y=log(ys);
c=log(cs);
i=log(is);
a=log(as);
r=log(rs);
R=Rs;
w=log(ws);
n=log(ns);
yn = log(ys/ns);
end;

shocks;
var ea = sda^2;
end;
resid(1);
steady;
check;

%Uses HP filter before computing moments.
stoch_simul(order =1, hp_filter =1600,periods=1000);
