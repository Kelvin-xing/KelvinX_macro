%NK with sticky price, Taylor rule
%we define the output gap as the difference btw output of sticky and
%flexible price 2014-04-11@ND
%we consider zero inflation steady state for a simple life here
%This file is written by Xiangyang Li;

var c i pi r n w mc a y vp pisharp x1 x2 yf outputgap m;
varexo ea ei;
parameters sigma beta psi eta phi epsilon theta phipi phiy rhoa rhoi sigmaa sigmam;
parameters cs is pis ns ws mcs as ys vps pisharps x1s x2s yfs outputgaps ms;

beta =.99;
sigma = 1;
eta = 1;
psi = 1;
epsilon =10;
theta = 1;
phipi = 1.5;
phiy = 0; %0.125/4;
rhoa = .95;
rhoi = 0.8;
sigmaa =.01;
sigmam=.01;
%the stickiness parameter
phi = .75;


%steady state calculation
is = 1/beta;
as = 1;
%zero inflation steady state
pis = 1;
pisharps= ((pis^(1-epsilon) - phi)/(1-phi))^(1/(1-epsilon));
vps = (1-phi)*(pis/pisharps)^epsilon/(1- pis^epsilon*phi);
mcs = (1-phi*beta*pis^epsilon)/(1-phi*beta*pis^(epsilon-1))
           *pis/pisharps*(epsilon-1)/epsilon;
ns = (vps^sigma*mcs/psi)^(1/(eta+sigma));
ys = as*ns/vps;
cs= ys;
yfs =((epsilon-1)/epsilon/psi)^(1/(sigma+eta))* as^((1+eta)/(sigma+eta));
outputgaps = ys/yfs;
ws=  as*mcs;
x1s = cs^(-sigma)*mcs*ys/(1-beta*phi*pis^epsilon);
x2s =cs^(-sigma)*ys/(1-beta*phi*pis^(epsilon-1));
ms = theta*is/(is -1)*cs^sigma;

model;
%(1) home Euler equation
exp(-sigma*c) = beta * exp(-sigma*c(+1))*exp(i)/exp(pi(+1));

%(2) labor supply
psi*exp(eta*n) =exp(-sigma*c) *exp(w);

%(3) labor demand
exp(mc) = exp(w)/exp(a);

%(4) Taylor interest rate rule
%actual output deviation
%i = (1-rhoi)*log(is) + rhoi*i(-1) + (1-rhoi)*(phipi*(pi - log(pis)) +phiy* (y - log(ys)) )+ ei;

% output gap
i = (1-rhoi)*log(is) + rhoi*i(-1) + (1-rhoi)*(phipi*(pi - log(pis)) 
         +phiy* (outputgap ))+ ei;

%(5) accounting identity
exp(c)=  exp(y);

%(6) the production technology
exp(y) = exp(a) *exp(n)/exp(vp);

%(7) the price dispersion
exp(vp) = (1-phi)*exp(-epsilon*pisharp)*exp(epsilon*pi) 
               + exp(epsilon*pi) *phi*exp(vp(-1));

%(8) inflation evolution
exp((1-epsilon)*pi) = (1-phi)*exp((1-epsilon)*pisharp) + phi;

%(9)the sticky price equation
exp(pisharp) = epsilon/(epsilon -1 )*exp(pi)*exp(x1)/exp(x2);

%(10) the auxiliary x1
exp(x1) = exp(-sigma*c)*exp(y)*exp(mc) +phi*beta*exp(epsilon*pi(+1))*exp(x1(+1));

%(11) the auxiliary x2
exp(x2) = exp(-sigma*c)*exp(y) +phi*beta*exp((epsilon-1)*pi(+1))*exp(x2(+1));

%(12) technology shock
a = rhoa*a(-1) + ea;

%(13)flexible output
exp(yf) = ((epsilon-1)/epsilon/psi)^(1/(sigma+eta))* exp(((1+eta)/(sigma+eta))*a);

%(14) output gap
outputgap = y - yf;

%(15) real interest rate
exp(r) = exp(i)/exp(pi(+1));

%(16) the real money balance
exp(m) = theta* exp(i)/(exp(i) -1 )*exp(sigma*c);

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
vp = log(vps);
pisharp = log(pisharps);
x1 = log(x1s);
x2 = log(x2s);
yf = log(yfs);
outputgap = log(outputgaps);
r =  log(is/pis);
m = log(ms);
end;

shocks;
var ea = .01^2;
var ei =.01^2;
end;

resid(1);
steady;
check;
%when used in loops, noprint and nograph is a good option.
%graph_format = none, only display,no save to disk
%since y and c is the same, we only plot y
stoch_simul(order=1) m i pi n w mc y yf r outputgap a;
