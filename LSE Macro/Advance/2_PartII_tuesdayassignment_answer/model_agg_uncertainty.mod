%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************

var c, k, z, Ka, r, w;
varexo e1,e2;

parameters k_ss,beta,nu,delta,zeta0,zeta1,zeta2,alpha,rho,b_0,b_K,b_z,sig_e1,sig_e2;

load parametervalues;
load aggregatelaw;
set_param_value('alpha',alpha)
set_param_value('nu',nu)
set_param_value('delta',delta)
set_param_value('zeta0',zeta0)
set_param_value('beta',beta)
set_param_value('rho',rho)
set_param_value('zeta1',zeta1)
set_param_value('zeta2',zeta2)
set_param_value('k_ss',k_ss)
set_param_value('b_0',b_0)
set_param_value('b_K',b_K)
set_param_value('b_z',b_z)
set_param_value('sig_e1',sig_e1)
set_param_value('sig_e2',sig_e2)

model;
Ka = b_0 + b_K*Ka(-1) + b_z*(z-1);

// rental rate capital
r = alpha*z*Ka(-1)^(alpha-1);

// wage rate
w = (1-alpha)*z*Ka(-1)^alpha;

// first-order condition for capital
c^(-nu)+zeta2-zeta1*exp(-zeta0*k)=beta*c(+1)^(-nu)*(r(+1)+1-delta);

// budget constraint
c+k=r*k(-1)+w*(1+e1)+(1-delta)*k(-1);

// law of motion aggregate productivity
z = (1-rho)+rho*z(-1)+e2;
end;

initval;
k = k_ss;
Ka = b_0/(1-b_K);
r = alpha*Ka^(alpha-1);
w = (1-alpha)*Ka^alpha;
z = 1;
c = w-delta*k;
end;
steady;

shocks;
var e1; stderr sig_e1;
var e2; stderr sig_e2;
end;

stoch_simul(order=2,nocorr,nomoments,IRF=0) k, Ka;
