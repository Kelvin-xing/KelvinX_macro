% This file is modified by Li Xiangyang @ 2012-5-6 to run in Dynare 4.2.0 or above
var c k lab z;
varexo e;

parameters bet the del alp tau rho s;

bet     = 0.987;
the     = 0.357;
del     = 0.012;
alp     = 0.4;
tau     = 2;
rho     = 0.95;
s       = 0.007;

model; 
% define model-local variables to improve readability;
# mar_c= (c^the*(1-lab)^(1-the))^(1-tau);
# mar_c1=(c(+1)^the*(1-lab(+1))^(1-the))^(1-tau);

%(1) Euler equation
mar_c/c=bet*(mar_c1/c(+1))*(1+alp*exp(z)*k(-1)^(alp-1)*lab^(1-alp)-del);

%(2) wage equation
c=the/(1-the)*(1-alp)*exp(z)*k(-1)^alp*lab^(-alp)*(1-lab);

%(3)capital accumulation equation
k=exp(z)*k(-1)^alp*lab^(1-alp)-c+(1-del)*k(-1);

%(4) technology shock
z=rho*z(-1)+e;
end;

initval;
k   = 1;
c   = 1;
lab = 0.3;
z   = 0;
e   = 0;
end;

shocks;
var e=s^2;
end;

steady;

%if periods not specify, there will be no simulations.
stoch_simul(periods=1000,irf=40,order=1); 

%save the simulated data to file
dynasave('simudata.mat');
