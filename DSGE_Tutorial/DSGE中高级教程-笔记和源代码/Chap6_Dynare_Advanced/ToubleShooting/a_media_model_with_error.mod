%This is the mod file 

var k kg i ig c lambda h w mc pi pistar R q Rk y H F mu x b; // 20
var tauc tauh tauk A G ; //5

varexo e_G, e_A, e_R,e_c,e_h,e_k e_ig e_mu e_x;

parameters ks kgs is igs cs lambdas hs ws mcs pis pistars Rs qs Rks ys Hs Fs mus xs As taucs tauhs tauks Gs bs;

parameters alpha beta delta gamma phi theta  thetab eta kappa     //Structure parameters
    rhoG psiGy psiGb        //fiscal rule
    rohig psiigy psiigb        //fiscal rule
    psicb psihb psikb
    psiRpi  psiRy
    rhomu rhoA rhoc rhoh rhok rhoig rhox;        //shocks autocorrelation

alpha=0.55;
beta=0.99;
delta=0.03;
gamma=5;
phi=0.75; 
theta=10;
eta=0.85;            
kappa=3.8;
thetab=0.19;
ig_y=0.03;

//fiscal policy efficients
rhoG=0.7;
psiGy=-0.4;
psiGb=-0.1;
rhoig=0.5;
psiigy=-0.2;                                        
psiigb=-0.095; 
psicb = 0;
psihb = 0;
psikb = 0;

//monetary policy coefficients
rhoR=0.75;
psiRpi=0.65;
psiRy=0.15;  

rhomu=0.35;
rhoA=0.3; 
rhoc=0.4; 
rhoh=0.6 ;
rhok=0.65;
rhox = 0.65;

//rhozeta=0.8; //param value of quantity rule and price rule following zhang(2009)
//psizetapi=-0.90; 
//psizetay=-0.5;


%steady states calculation
pis=1;    
pistars =1;
As = 1;
qs =1;   
mus = 1;
xs = 1;

taucs=0.17;
tauhs=0.2;
tauks=0.2; 

Rs=pis/beta;
mcs=(theta-1)/theta;
Rks=(Rs-1+delta)/(1-tauks);

c_y=0.45;                               
lambday=1/((1-eta)*(1+taucs)*c_y);       
ig_kg=delta;
kg_y=ig_y/ig_kg; 
ws=(alpha^(alpha/(1-alpha-thetab)))*((1-alpha)^(1+thetab/(1+gamma-alpha-thetab-gamma*alpha-gamma*thetab)))*((lambday*(1-tauhs))^(thetab/(1+gamma-alpha-thetab-gamma*alpha-gamma*thetab)))*((kg_y)^(thetab/(1-alpha-thetab)))*(mcs^((1+gamma-gamma*thetab)/(1+gamma-alpha-thetab-gamma*alpha-gamma*thetab)))/(Rks^(alpha/(1-alpha-thetab)));

h_y=(1-alpha)*mcs/ws;
k_y=alpha*mcs/Rks;
i_k=delta;
i_y=i_k*k_y;
PI_y=1-mcs;

m_y=1/(lambday*(1-beta));                         
G_y=1-c_y-i_y-ig_y;                     
    
//error usage,Rk should be Rks
               
b_y=(c_y*taucs+ws*h_y*tauhs+Rk*k_y*tauks-G_y-ig_y)/(Rs-1); 
ys=((lambday*(1-tauhs))^(1/(1+gamma)))*ws/(((1-alpha)*mcs)^(gamma/(1+gamma)));       
ks=k_y*ys;       
cs = c_y*ys;
is = i_y*ys;
Gs =G_y*ys;
kgs = kg_y*ys;
igs = ig_kg*kgs;

lambdas = hs^gamma/ws/(1-tauhs);
hs = h_y*ys;
Hs = lambdas*ys*mcs/(1-phi*beta);
Fs= lambdas*ys/(1-phi*beta);
bs = (cs*taucs + ws*hs*tauhs +Rks*ks*tauks-Gs - igs )/(Rs/pis -1);

model;
#inv = 1- kappa/2*(1- exp(i - i(-1)))^2 +kappa*(1- exp(i - i(-1)))*exp(i - i(-1));
#inv1 = (1- exp(i(1) - i))*exp(i(1) - i)^2;
//(1) captial stock evolution
exp(k)=(1-delta)*exp(k)+exp(x)*exp(i)*(1-kappa/2*(exp(i-i(-1))-1)^2);

//(2) foc of consumption
exp(mu)/(exp(c)-eta*exp(c(-1))) = exp(lambda)*(1+exp(tauc));

//(3) foc of labor
exp(gamma*h) = exp(lambda+w)*(1-exp(tauh));

//(4) foc of bonds
exp(lambda) = beta*exp(lambda(1)-pi(1)+R);

//(5) foc of capital
exp(lambda+q) = beta*exp(lambda(1))*(exp(Rk(1))*(1-exp(tauk(1)))+exp(q(1))*(1-delta));

//(6) foc of investment
exp(lambda)*(1-exp(q+x)*inv)-beta*kappa*exp(lambda(1)+q(1)+x(1))*inv1;

//(7) production fucntion
exp(y) = exp(A+alpha*k(-1)+(1-alpha)*h+thetab*kg(-1));

//(8) government capital
exp(kg) = (1-delta)*exp(kg(-1)) + exp(ig);

//(9) labor demand
exp(h) = (1-alpha)*exp(y+mc -w);

//(10)capital demand
exp(k) = alpha*exp(y+mc -Rk);

//(11) resource constraint
exp(y) = exp(c) + exp(i) + exp(ig) + exp(G);

//(12) government budget constraint
exp(G) + exp(ig) + exp(R(-1) + b(-1) - pi) = exp(tauc + c) + exp(w+h + tauh) + exp(Rk+k(-1)+tauk) + exp(b);

//(13) the price evolution 
exp((1-theta)*pi) = (1-phi)*exp((1-theta)*pistar) + phi;

//(14) the optimal price
exp(pistar) = theta/(theta -1 )*exp(pi)*exp(H)/exp(F);

%(15) the auxiliary H
exp(H) = exp(lambda)*exp(y)*exp(mc) +phi*beta*exp(theta*pi(+1))*exp(H(+1));

%(16) the auxiliary F
exp(F) = exp(lambda)*exp(y) +phi*beta*exp((theta-1)*pi(+1))*exp(F(+1));

%(17) technology shock
A = rhoA*A(-1) + e_A;

//(18) fiscal rule 1 
// error using level and log-level together;
G - Gs = rhoG*(G(-1) - Gs) +psiGy*(y(-1) - ys) + psiGb*(b(-1) -bs) + e_G;

//(19) fiscal rule 2
ig - log(igs) = rhoig*(ig(-1) - log(igs)) + psiigy*(y(-1) - log(ys)) + psiigb*(b(-1) - log(bs)) + e_ig;

//(20) monetary policy
R - log(Rs) = rhoR*(R(-1) - log(Rs)) +psiRpi*pi(1) + psiRy*(y - log(ys)) +e_R;

//(21) mu motion
mu - log(mus) = rhomu*(mu(-1) - log(mus)) +e_mu;

//(22) x motion
x - log(xs) = rhox*(x(-1) - log(xs)) +e_x;

//(23) tauc motion
tauc - log(taucs) = rhoc*(tauc(-1) - log(taucs)) +psicb*(b(-1) - log(bs))+e_c;

//(24) tauh motion
tauh - log(tauhs) = rhoh*(tauh(-1) - log(tauhs)) +psihb*(b(-1) -log(bs))+e_h;

//(25) tauk motion
tauk - log(tauks) = rhok*(tauk(-1) - log(tauks)) +psikb*(b(-1) -log(bs))+e_k;
end;

initval;
k = log(ks);
kg = log(kgs);
i = log(is);
ig = log(igs);
c = log(cs);
lambda = log(lambdas);
h = log(hs);
w = log(ws);
mc = log(mcs);
pi = log(pis);
pistar = log(pistars);
R = log(Rs);
q = log(qs);
Rk = log(Rks);
y = log(ys);
H = log(Hs);
F = log(Fs);
mu = log(mus);
x = log(xs);
b = log(bs);
tauc = log(taucs);
tauh = log(tauhs);
tauk = log(tauks);
A = log(As);
G = log(Gs);
end;

//syntax error
shock; 
var e_A = 0.01^2;
var e_G = 0.01^2;
var e_ig =0.01^2;
var e_R  = 0.01^2;
var e_c = 0;
var e_k = 0;
var e_h  =0;
var e_mu =0;
var e_x  =0;
end;

resid(1);
steady;
check;

stoch_simul(order=1);
