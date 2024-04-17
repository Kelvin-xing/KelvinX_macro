% we have 13 endogenous variables
var kbar //capital stock
 i           // investment
c           //consumption
R         // nominal rate
omegabar // cutoff value of risk
Rk          // return rate of capital
n            //net worth of enterpreneur
sigma   // variable standard devation of risk shock
spread  // interest rate spread
credit   //  loan from bank by enterpreneur
bankrupt  //rate of bankruptancy
GDP     // output
wedge   //financial wedge
;

varexo sigma_e;

parameters sigma_ss mu gamma alpha delta beta rhosigma gam;

sigma_ss=0.26; %steady state of standard dev. of risk shock
mu=0.21;          %monitoring cost parameter
gamma=0.97;    %fracation of the income goes to the households
alpha=1/3;         %capital share
delta=0.02;        %depreciation rate
beta=1.03^(-.25); %discount rate
rhosigma=0.97; %persistence
gam=1;             %log utility

model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary expressions.  These simplify the equations without adding
% additional variables. Noticing the timing scprits of the variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # z        = (log(omegabar) + sigma(-1)^2 / 2) / sigma(-1);
  # zplus   = (log(omegabar(+1)) + sigma^2 / 2) / sigma;
  # F          = normcdf(z);
  # G          = normcdf(z - sigma(-1));
  # d          =  mu * G  * (1 + Rk) * kbar(-1);
  # Fp1        = normcdf(zplus);
  # Gp1        = normcdf(zplus - sigma);
  # GAMMA      = omegabar*(1-F)+G;
  # GAMMAp1    = omegabar(+1)*(1-Fp1)+Gp1;
  # dFp1       = normpdf(zplus) / omegabar(+1) / sigma;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1) resource constraint
    kbar(-1)^alpha = c + i + d ;
  
% (2)Law of motion for capital, capital producer technolgy
    kbar = (1 - delta) * kbar(-1) + i;
  
% (3) Household FOC w.r.t. risk-free bonds
    1 = beta * (c(+1)/c)^(-gam)*(1+R);
  
% (4) efficiency condition associated with the contract
    (1-Fp1)/(1-GAMMAp1)-((1 + Rk(+1)) / (1 + R))*(1-Fp1-mu*omegabar(+1)*dFp1)
    /(1-((1 + Rk(+1)) / (1 + R))*(GAMMAp1-mu*Gp1));
  
% (5) Definition of return of entrepreneurs, Rk
    1 + Rk = alpha*kbar(-1)^(alpha-1)  + 1 - delta;
    
%(6) Zero bank profit condition, i.e. the definition of leverage ratio
    kbar(-1)*(1 - ((1+Rk)/(1+R(-1)))*(GAMMA-mu*G)) - n(-1);
  
%(7) Law of motion of net worth, omitting the entrepreneurial wage
    n = gamma*(1-GAMMA)*(1+Rk)*kbar(-1);

%(8) borrowing:
    credit = kbar-n;

%(9) interest rate spread
    spread = ((1+Rk)/(1+R(-1)))*(kbar(-1)/credit(-1))*omegabar;

%(10) Shock 
   log(sigma / sigma_ss) = rhosigma * log(sigma(-1) / sigma_ss)  + sigma_e ;

%(11) the output of the economy
   GDP = c + i;

%(12) the financial wedge
   wedge = (1-(1+R)/(1+Rk(+1)));

%(13) bankruptcy rate
   bankrupt = F;
end;

shocks;
var sigma_e; 
stderr .11;
end;

%this directive will invoke steady state file to have s.s.
steady;

stoch_simul(periods=1000,order=1,irf=20) kbar, i, c, R, omegabar, 
                   Rk, n, sigma, spread, credit, bankrupt, GDP, wedge;
