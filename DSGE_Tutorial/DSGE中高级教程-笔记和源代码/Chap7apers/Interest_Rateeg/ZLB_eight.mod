%Zero lower bound, interest rate peg,NK with sticky price, Taylor rule
%2014-04-11@ND
%we consider zero inflation steady state for a simple life here
%This file is written by Xiangyang Li, xli21@nd.edu;

var c i pi r n w mc a y vp yf pisharp x1 x2 g nex s1 s2 s3 s4 s5 s6 s7 s8;
varexo ea ei eg ex;
parameters sigma beta psi eta phi epsilon phipi phiy etag etax;
parameters rhog rhox rhoa rhoi sigmaa sigmai sigmax sigmag;
parameters cs is pis ns ws mcs as ys yfs vps pisharps x1s x2s gs nexs;

beta =.99;
sigma = 1;
eta = 1;
psi = 1;
epsilon =10;
phipi = 1.5;
phiy = 0;
etag = .1413; % 1993-2012,annual data,CEIC
etax = .0320; % 2005-2012 annual data,0.0459,CEIC;1994-2013,0.0320
rhoa = .95;
rhoi = .0;
rhog =  .9658;
rhox =  .8257;

sigmaa =.1;%shut off
sigmai=.1; %shut off 
sigmax =.5010; 
sigmag = .1459; 


%the stickiness parameter
phi = .75;
%zero inflation steady state
pis = 1;

%steady state calculation
is = 1/beta*pis;
as = 1; %the steady state value of technology variable

pisharps= ((pis^(1-epsilon) - phi)/(1-phi))^(1/(1-epsilon));
vps = (1-phi)*(pis/pisharps)^epsilon/(1- pis^epsilon*phi);
mcs = (1-phi*beta*pis^epsilon)/(1-phi*beta*pis^(epsilon-1))*pisharps/pis*(epsilon-1)/epsilon;
ns = (vps^sigma*mcs/psi/(1-etax-etag)^(sigma))^(1/(eta+sigma));
ys = as*ns/vps;
yfs =  ((epsilon-1)/epsilon/psi/(1-etag - etax)^sigma)^(1/(sigma + eta))*as^((1+eta)/(sigma + eta));
gs = ys*etag;
nexs = ys*etax;
cs= ys*(1-etag - etax);
ws=  as*mcs;
x1s = cs^(-sigma)*mcs*ys/(1-beta*phi*pis^epsilon);
x2s =cs^(-sigma)*ys/(1-beta*phi*pis^(epsilon-1));

model;
%(1) home Euler equation
exp(-sigma*c) = beta * exp(-sigma*c(+1))*exp(i)/exp(pi(+1));

%(2) labor supply
psi*exp(eta*n) =exp(-sigma*c) *exp(w);

%(3) labor demand
exp(mc) = exp(w)/exp(a);

%(4) accounting identity
exp(c) +exp(g) + exp(nex) =  exp(y);

%(5) the production technology
exp(y) = exp(a) *exp(n)/exp(vp);

%(6) the price dispersion
exp(vp) = (1-phi)*exp(-epsilon*pisharp)*exp(epsilon*pi) + exp(epsilon*pi) *phi*exp(vp(-1));

%(7) inflation evolution
exp((1-epsilon)*pi) = (1-phi)*exp((1-epsilon)*pisharp) + phi;

%(8)the sticky price equation
exp(pisharp) = epsilon/(epsilon -1 )*exp(pi)*exp(x1)/exp(x2);

%(9) the auxiliary x1
exp(x1) = exp(-sigma*c)*exp(y)*exp(mc) +phi*beta*exp(epsilon*pi(+1))*exp(x1(+1));

%(10) the auxiliary x2
exp(x2) = exp(-sigma*c)*exp(y) +phi*beta*exp((epsilon-1)*pi(+1))*exp(x2(+1));

%(11) technology shock
a = rhoa*a(-1) + ea;

%(12)government spending shock
g = (1-rhog)*log(gs) + rhog*g(-1) + eg;

%(13) net export shock
nex = (1-rhox)*log(nexs) +  nex(-1)*rhox + ex;

%(14) real interest rate
exp(r) = exp(i)/exp(pi(+1));

%(15) the flexible price output
exp(yf) = ((epsilon-1)/epsilon/psi/(1-etag - etax)^sigma)^(1/(sigma + eta))*exp(a*(1+eta)/(sigma+eta));

%(16 -24) for interest rate peg model,H=8 
i = s1(-1);
s1 = s2(-1);
s2 = s3(-1);
s3 = s4(-1);
s4 =s5(-1);
s5 =s6(-1);
s6 =s7(-1);
s7 = s8(-1);
s8 =  (1-rhoi)*log(is) + rhoi*s8(-1) + (1-rhoi)*(phipi*(pi(+8) - log(pis)) + phiy*(y(+8)  - yf(+8))) +ei(+8); 

end;
initval;
c = log(cs);
i = log(is);
pi = log(pis);
n = log(ns);
w = log(ws);
mc = log(mcs);
a = log(as);
y = log(ys);
yf = log(yfs);
vp = log(vps);
pisharp = log(pisharps);
x1 = log(x1s);
x2 = log(x2s);
g =  log(gs);
r =  log(is/pis);
nex = log(nexs);
s1 = log(is);
s2 = log(is);
s3 = log(is);
s4 = log(is);
s5 = log(is);
s6 = log(is);
s7 = log(is);
s8 = log(is);
end;

shocks;
var ea = sigmaa^2;
var ei = sigmai^2;
var eg = sigmag^2;
var ex = sigmax^2;
end;

%resid(1);
steady;
check;
%when used in loops, noprint and nograph is a good option.
%graph_format = none, only display,no save to disk
%since y and c is the same, we only plot y
stoch_simul(order=1,nograph) i pi n w mc y r g c nex yf;
