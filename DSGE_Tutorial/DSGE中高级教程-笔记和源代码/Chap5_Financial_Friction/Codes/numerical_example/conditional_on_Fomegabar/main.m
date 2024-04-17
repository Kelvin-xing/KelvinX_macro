clear all
close all
dbstop if error

R=1.011470654984364;
Rk=1.018806396097570;
sp=Rk/R;%risk spread
mu=0.2148945111281970;%monitor cost
Fomegabar=0.0055885692972290; %bankruptancy rate are set.

%when given Fomega,mu and spread, solve omega, then sigma
[sigma,omega,Gamma,Gam_muG,ixx] = get_omega_cond_Fomegabar(sp,Fomegabar,mu);


leverage=1/(1-Gam_muG*sp); %bank zero profit condition
spread=sp*omega*leverage/(leverage-1);%interest rate  spread, gross
spreadAR=400*(spread-1); %annulized interest rate spread;
util=(1-Gamma)*sp*leverage;


