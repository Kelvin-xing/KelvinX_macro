function [ff,sp,omega,G,F,Gamma,Gam_muG,Fprime,k,n]=get_bank_zero_profit_diff_cond_on_Rk(Rk,R,sigma,mu,alpha,gamma,delta,z)
%starting from efficiency condition, when given spread(sp1,since R is known then Rk is known), then finding a
%omega from efficiency condition. then find  k and n. Finally, finding
%differnce of the two sides of the bank zero conditions.
sp=Rk/R;
%when given sigma, mu and spread(three assumptions, i.e.,three calibration), 
%'get_omega_cond_on_sigma.m' finds an omega that solves the efficiency
%condition, i.e., the foc of the entrepreneurial utility maximization
%problem.
[omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_on_sigma(sigma,mu,sp);

%rate of return on the capital
k=(z*alpha/(Rk-(1-delta)))^(1/(1-alpha));

%aggregate net worth, and  w_e is ignored
n=gamma*(1-Gamma)*Rk*k;

%bank zero profit conditions
ff=(k/n)-1/(1-sp*Gam_muG);