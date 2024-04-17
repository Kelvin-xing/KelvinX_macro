//Written by Xiangyang@BJ,2015-10
//the log-level model of the CGG;
var c k a f;
varexo eps_a;

parameters gamma alpha delta beta rho Veps;
parameters cs ks as fs;

beta = 0.99;
gamma = 2; //20
alpha = 0.36;
delta =0.02;
rho =0.95;
Veps = 0.01^2;
ks = (alpha*beta/(1-(1-delta)*beta))^(1/(1-alpha));
as = 0;
fs = ks^alpha+(1-delta)*ks;
cs  = fs - ks;

model;
//(1) resource constraint
exp(c) + exp(k) = exp(f);

//(2) the production technology
exp(f) = exp(a)*exp(alpha*k(-1)) + (1-delta)*exp(k(-1));

//(3) the Euler equation
beta*exp(-gamma*c(+1))*(alpha*exp((alpha-1)*k)
*exp(a(+1))+(1-delta))=exp(-gamma*c);

//(4) the technology motion
a = rho*a(-1) + eps_a;
end;

initval;
a = as;
f = log(fs);
k = log(ks);
c = log(cs);
end;

shocks;
var eps_a = Veps;
end;

//solve and simulate the model, set qz_zero_threshold=1e-15 if gamma=20
stoch_simul(order =2,nograph,qz_zero_threshold=1e-15);


