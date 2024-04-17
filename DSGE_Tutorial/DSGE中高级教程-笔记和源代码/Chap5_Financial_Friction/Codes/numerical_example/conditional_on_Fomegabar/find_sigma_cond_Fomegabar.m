function [sigma] = find_sigma_cond_Fomegabar(omega,Fomega)
% Given omega and Fomega, This routine finds the value of sigma for the lognormal distribution
% with mean forced to equal unity, such that prob < omega = Fomega.
%For Fomega = .03, this program requires
%omega to live in the interval [0.000000001,0.999];
sigma1=.00000001;
ff1 = logncdf(omega,-sigma1^2/2,sigma1);
sigma2=5;
ff2 = logncdf(omega,-sigma2^2/2,sigma2);

if ff1*ff2 > 0
    error('fatal (ffsigma) failed to bracket sigma')
end

%give a relatively large interval for sigma;
x0=[sigma1 sigma2];
opt=optimset('diagnostics','on','TolX',1e-12,'TolFun',1e-12);
[sigma,fval,exitflag,output] =   fzero(@find_lognormal_difference,x0,opt, ...
    omega,Fomega);
if abs(fval ) > .1e-9 || abs(imag(sigma))>.1e-10 || sigma < 0 || exitflag <= 0
    error('fatal (ffsigma) failed to find sigma')
end

    
    

