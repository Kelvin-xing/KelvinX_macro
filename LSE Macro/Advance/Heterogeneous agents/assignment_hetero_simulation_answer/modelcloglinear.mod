%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
% modelcloglinear.mod
% standard growth model
% variables are the log of consumption, the log of capital, and the log of productivity
% Dynare generates a law of motion that is linear in these variables (since order = 1)
 
var lnc, lnk, lnz;
varexo e;

parameters beta, rho, alpha, nu, delta, sig;

alpha = 0.2;
rho   = 0.95;
beta  = 0.99;
nu    = 1;
sig   = 0.007;
delta = 0.10;
model;
exp(-nu*lnc)=beta*(exp(-nu*lnc(+1)))*(exp(lnz(+1))*alpha*exp((alpha-1)*lnk)+1-delta);
exp(lnc)+exp(lnk)=exp(lnz+alpha*lnk(-1))+(1-delta)*exp(lnk(-1));
lnz = rho*lnz(-1)+e;
end;

initval;
lnk = (1/(1-alpha))*log(beta*alpha/(1-beta*(1-delta)));
lnc = log(exp(alpha*lnk)-delta*exp(lnk));
lnz = 0;
end;

shocks;
var e; stderr sig;
end;

stoch_simul(order=1,nocorr,nomoments,IRF=0) lnk;