%CIA 
var c n w R i k a y I dlnm pi m r mu lambda;
varexo ea em;
parameters theta beta alpha delta rhom pistar sigmam sigmaa rhoa;
parameters cs ns ws Rs is ks ys Is ms rs mus lambdas;
beta =.99;
alpha =1/3;
delta=.025;
rhom=.5;
rhoa = .5;
sigmam =.01;
sigmaa=.01;
pistar = 1.02;
Rs = 1/beta - 1+delta;
kn = (alpha/Rs)^(1/(1-alpha));
ws = (1-alpha)*kn^alpha;
is = pistar/beta;
ns = 1/3;
ks = kn*ns;
Is = delta*ks;
ys = kn^alpha*ns;
cs = ys - Is;
ms = cs*pistar;
lambdas = beta/cs/pistar;
mus = 1/cs - lambdas;
rs = is/pistar;
theta = ws*lambdas*(1-ns);

model;
%(1) labor supply equation
theta/(1-exp(n)) = exp(lambda)*exp(w);

%(2) Euler equation
exp(lambda) = beta*(exp(lambda(+1))* ( exp(R(+1)) + 1-delta) );

%(3) bonds 
exp(lambda) = beta*(exp(lambda(+1))* ( exp(i)/exp(pi(+1)) ));

%(4) capital returns
exp(R) = alpha*exp(a)*exp(k(-1))^(alpha -1)*exp(n)^(1-alpha);

%(5) wage
exp(w) = (1-alpha)*exp(a)*exp(k(-1))^alpha*exp(n)^(-alpha);

%(6) production technology
exp(y) = exp(a)*exp(k(-1))^alpha*exp(n)^(1-alpha);

%(7) capital accumulation
exp(k) = exp(I) + exp(k(-1))*(1 - delta );

%(8) GDP identity
exp(y) = exp(I) + exp(c);

%(9) Fisher relationship
exp(r) = exp(i)/exp(pi(+1));

%(10) real money balance
exp(lambda) = beta*(exp(mu(+1)) +exp(lambda(+1)))/exp(pi(+1));

%(11) monetary shock
dlnm = (1-rhom)*log(pistar) - (pi) +rhom*(pi(-1)) + rhom*dlnm(-1) +em;

%(12) technology shock
a = rhoa*a(-1) + ea;

%(13) price level and money level
dlnm = m - m(-1);

%(14) marginal consumption
1/exp(c) = exp(lambda) + exp(mu);

%(15) CIA
exp(m(-1)) = exp(c)*exp(pi);

end;

initval;
c = log(cs);
n = log(ns);
w = log(ws);
R = log(Rs);
i = log(is);
k = log(ks);
a = 0;
y = log(ys);
I = log(Is);
dlnm = 0;
pi = log(pistar);
m = log(ms);
r = log(rs);
lambda = log(lambdas);
mu = log(mus);
end;

resid;
steady;

shocks;
var ea = sigmaa^2;
var em = sigmam^2;
end;
stoch_simul(order=1);
