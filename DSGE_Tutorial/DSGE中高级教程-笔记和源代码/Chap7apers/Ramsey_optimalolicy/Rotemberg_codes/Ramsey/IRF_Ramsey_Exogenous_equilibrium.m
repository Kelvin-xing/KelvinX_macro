%this file will plot the RBC model's IRF from exogenous monetary policy
%rule and Ramsey equilibrium.

clear all;
close all;
clc;

%run the dynare mod file and store the results;
dynare exogenous_monetary_rule noclearall;
save exogenous_mon_rule oo_ M_ options_;

dynare Rotemberg_ramsey_policy noclearall;
save ramsey_policy oo_ M_ options_

%plot the IRF together

