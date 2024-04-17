%NK with sticky price 
%we define the output gap as the difference btw output of sticky and
%flexible price 2014-04-11@ND
%we consider zero inflation steady state for a simple life here
%This file is written by Xiangyang Li;
%The defintion of nominal interest rate is slightly different from the note.

var c i pi n w m mc a y vp pisharp x1 x2 dlnm yf outputgap r;
varexo ea em;
parameters sigma beta psi eta phi epsilon gamma rhoa rhom sigmaa sigmam;
parameters cs is pis ns ws ms mcs as ys vps pisharps x1s x2s dlnms yfs outputgaps rs;

beta =.99;
sigma = 1;
eta = 1;
psi = 1;
epsilon =10;
gamma = 1;
rhoa = .95;
rhom = 0;
sigmaa =.01;
sigmam=.01;
%the stickiness parameter
phi = .75;

%steady state calculation
is = 1/beta;
rs = is - 1; %real interest rate;
as = 1;
dlnms = 0;
%zero inflation steady state
pis = 1;
pisharps= ((pis^(1-epsilon) - phi)/(1-phi))^(1/(1-epsilon));
vps = (1-phi)*(pis/pisharps)^epsilon/(1- pis^epsilon*phi);
mcs = (1-phi*beta*pis^epsilon)/(1-phi*beta*pis^(epsilon-1))
           *pis/pisharps*(epsilon-1)/epsilon;
ns = (vps^sigma*mcs/psi)^(1/(eta+sigma));
ys = as*ns/vps;
cs= ys;
ms = gamma*is/(is-1)*ys^sigma;
yfs =((epsilon-1)/epsilon/psi)^(1/(sigma+eta))* as^((1+eta)/(sigma+eta));
outputgaps = ys/yfs;
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

%(4) money demand
exp(m) = gamma*exp(i)*exp(sigma*c)/(exp(i) -1 );

%(5) money supply
dlnm = (1-rhom)*(pis-1) - pi + rhom*dlnm(-1) + rhom*pi(-1) + em;

%(6) real money balance growth 
dlnm = m - m(-1);

%(7) accounting identity
exp(c)=  exp(y);

%(8) the production technology
exp(y) = exp(a) *exp(n)/exp(vp);

%(9) the price dispersion
exp(vp) = (1-phi)*exp(-epsilon*pisharp)*exp(epsilon*pi) 
                + exp(epsilon*pi) *phi*exp(vp(-1));

%(10) inflation evolution
exp((1-epsilon)*pi) = (1-phi)*exp((1-epsilon)*pisharp) + phi;

%(11)the sticky price equation
exp(pisharp) = epsilon/(epsilon -1 )*exp(pi)*exp(x1)/exp(x2);

%(12) the auxiliary x1
exp(x1) = exp(-sigma*c)*exp(y)*exp(mc) 
                +phi*beta*exp(epsilon*pi(+1))*exp(x1(+1));

%(13) the auxiliary x2
exp(x2) = exp(-sigma*c)*exp(y) +phi*beta*exp((epsilon-1)*pi(+1))*exp(x2(+1));

%(14) technology shock
a = rhoa*a(-1) + ea;

%(15)flexible output
exp(yf) = ((epsilon-1)/epsilon/psi)^(1/(sigma+eta))* exp(((1+eta)/(sigma+eta))*a);

%(16) output gap
outputgap = y - yf;

%(17) real rate eq., Fisher equation
r = exp(i) - exp(pi(+1));
end;
initval;
c = log(cs);
i = log(is);
pi = log(pis);
n = log(ns);
w = log(ws);
m = log(ms);
mc = log(mcs);
a = log(as);
y = log(ys);
vp = log(vps);
pisharp = log(pisharps);
x1 = log(x1s);
x2 = log(x2s);
dlnm = dlnms;
yf = log(yfs);
outputgap = log(outputgaps);
r = is - pis;
end;

shocks;
var ea = .01^2;
var em =.01^2;
end;

resid(1);
steady;
check;
%when used in loops, noprint is a good option.
stoch_simul(order=1) i pi n r m mc a y yf outputgap;
