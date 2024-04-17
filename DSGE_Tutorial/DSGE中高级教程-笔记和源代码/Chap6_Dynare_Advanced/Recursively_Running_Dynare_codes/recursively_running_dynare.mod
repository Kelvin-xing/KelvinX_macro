%This program is written by Xiangyang Li @UND
%This is a very simple RBC model;

var y i r n w k a R c;
varexo e; 
parameters alpha beta delta rho sigma psi eta sigmae;
parameters yss iss rss nss wss kss Rss css;

%load the parameters saved in the mat file;
%load is the matlab command which can be used here;
load parametersaved;

%this function will pass the saved parameter into this mod file
%set_param_value() is the built-in function of dynare.
set_param_value('alpha',alpha);

beta = .99;
delta = .02;
rho = .97; 
sigma = 1; 
psi =3; 
eta = 1;
sigmae = .01;

rss = 1/beta - 1;
Rss= 1/beta - 1 + delta;
kn = (alpha/Rss)^(1/(1-alpha));
wss = (1-alpha)*kn^alpha;
nss = ((kn^alpha - delta*kn)*(psi/wss)^(1/sigma))^(-1/(1+eta/sigma));
kss = kn*nss;
yss = kn^alpha*nss;
iss=  delta*kss;
css = yss - iss;

model;
% (1) Euler equation, capital
exp(c)^(-sigma)=beta*exp(c(+1))^(-sigma)*(exp(R(+1))+(1-delta));

% (2) Euler equation, bonds
exp(c)^(-sigma)=beta*(1+r)*exp(c(+1))^(-sigma);

% (3) Labor supply
psi*exp(n)^(eta)=exp(c)^(-sigma)*exp(w);

% (4) Production func
exp(y)=exp(a)*exp(k(-1))^(alpha)*exp(n)^(1-alpha);

% (5) Capital demand
exp(R)=alpha*exp(a)*exp(k(-1))^(alpha-1)*exp(n)^(1-alpha);

% (6) Labor demand
exp(w)=(1-alpha)*exp(a)*exp(k(-1))^(alpha)*exp(n)^(-alpha);

% (7) Resource constraint
exp(y)=exp(c)+exp(i);

% (8) Capital accumulation
exp(k)=exp(i)+(1-delta)*exp(k(-1));

% (9) Productivity shock
a=rho*a(-1)+e;
end;

initval;
k=log(kss);
y=log(yss);
c=log(css);
i=log(iss);
a=0;
r=rss;
R=Rss;
w=log(wss);
n=log(nss);
end;

shocks;
var e = sigmae^2;
end;
resid(1); 
steady;
check;

%everything will show up in each loop which will consume a lot of time
%stoch_simul(order =1);

%useful for loops, most info. will not display,
%but compling info. still show up.
%stoch_simul(order =1,nograph,nofunctions);
stoch_simul(order =1,noprint,nograph);