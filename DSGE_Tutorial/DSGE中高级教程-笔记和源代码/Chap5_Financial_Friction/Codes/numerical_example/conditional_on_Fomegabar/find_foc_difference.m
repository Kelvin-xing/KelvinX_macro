function [ff,ix,Gamma,Gam_muG,sigma] = find_foc_difference(omega,risk_spread,Fomegabar,mu)

% Given omega and Fomeag, this function compute the Gamma, Gamma - mu*G 
% and the difference between the two sides of FOC. 
% An indicator will also return to indicate if there
% are something at odd with the foc in utility maximization problem.

if risk_spread<=1,error('(findomega) spread must be larger than 1'),end

%find a sigma such that logncdf(omega,-sigma^2/2,sigma)==Fomegabar when
%omega given.
sigma   =   find_sigma_cond_Fomegabar(omega,Fomegabar);

%for given sigma and omega,mu,Fomegabar, to find the difference between the
%two sides of the FOC.
% The following definitions are taken from BGG page 52 (NBER WP 6455)
%simplifying parameter
z       =   (log(omega)+sigma^2/2)/sigma; 
%Gamma(omega), the share of entrepreneurial profits 
%given to bank with the monitoring costs
Gamma   =   normcdf(z-sigma)+omega*(1-normcdf(z)); 

%Gamma - mu*G, share of entrepreneurial profits given 
%to bank net of monitoring costs
Gam_muG =   (1-mu)*normcdf(z-sigma)+omega*(1-normcdf(z)); 

% ff stands for the difference of the two sides of the FOC,i.e. the foc of entrepreneurial utility maximization, when
% ff==0, then omegabar is found.
ff  =   (1-Gamma)*risk_spread+(1-Fomegabar)/(1-Fomegabar-mu*omega*lognpdf(omega,-sigma^2/2,sigma))*(risk_spread*Gam_muG-1);
%an indicator, if equal to zero, nothing wrong. otherwise, something
%happened
ix=0;

%if ff complex number or parameters passed in are odd with FOC of utility
%maximization problem, see note for more information.
%indicator will turn on at this time.
if max(abs(imag(ff))) > .1e-8 || (1-Fomegabar)/(1-Fomegabar-mu*omega*lognpdf(omega,-sigma^2/2,sigma)) < 0
    ix=1;
end