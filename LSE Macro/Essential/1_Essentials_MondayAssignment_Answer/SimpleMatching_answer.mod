%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part I: The Essentials
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%==========================================================================
%                       Simple matching model
%==========================================================================

var c, n, v, y, m, w, Q, g, pf, z; 

varexo e;

parameters beta, nu, phi, psi, mu, rho, rhox, sigma, omega, omega_e, pf_ss, n_ss, g_ss;

beta    = 0.96;
rho     = 0.98;
nu      = 1;
rhox    = 0.027;
phi     = 0.3;
mu      = 0.5;
sigma   = 0.007;
omega   = 0.8;      %If omega = 1, wages are flexible; if omega = 0, wages are fixed

% loading entrepreneur share parameter
load ent_share_value;
set_param_value('omega_e',om_e);

%Calibrated endogenous variable:
pf_ss = 0.338;

%Solving for a parameter psi:
n_ss  = 1/((pf_ss/phi)^(1/mu)*(rhox/pf_ss)+1);
g_ss  = beta*(1/(1-beta*(1-rhox)))*(omega_e);
psi   = pf_ss*g_ss;

model;
%Household bc:
c  + psi*v = w*n(-1) + Q;

%Firm:
Q   = y - w*n(-1);
y   = exp(z)*n(-1);
n   = (1-rhox)*n(-1) + pf*v;
g   = beta*(c(+1)/c)^(-nu)*(exp(z(+1)) - w(+1) + (1-rhox)*g(+1));
psi = pf*g;

%Sharing rule and wage:
w  = (1-omega_e)*(omega*exp(z) + (1-omega));

%Matching:
m   = phi*(1-n(-1))^mu*v^(1-mu);
pf  = m/v;

%Shock process:
z   = rho*z(-1) + e;

end;

initval;
pf = pf_ss;
n  = n_ss;
g  = g_ss;
y  = n;
v  = rhox*n/pf;
m  = pf*v;
w  = 1-omega_e;
Q  = n - w*n;
c  = w*n + Q - psi*v;
z = 0;
 
end;

steady;

check;

shocks;
var e; stderr sigma;
end;

stoch_simul(order=1,nocorr,nomoments,nograph,periods=10000) c n v y m w Q g pf z;