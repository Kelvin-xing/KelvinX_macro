function [ff,G,F,Gamma,Gam_muG,Fprime]=find_foc_difference(omega,sig,mu,sp)

% this function only serves to calculate some expression in the foc
% and return the difference of the two sides of foc from utility maximization.
z       =   (log(omega)+sig^2/2)/sig;

%the function G and F defined in the note
F=logncdf(omega,-sig^2/2,sig);
G      =   normcdf(z-sig);

Gamma   =   G+omega*(1-F);
Gam_muG =   (1-mu)*G+omega*(1-F);
Fprime  =   (1/sqrt(2*3.14159265358979))*exp(-z^2/2)/(omega*sig);
%Fprime can also be defined as the following
%Fprime=lognpdf(omega,-sig^2/2,sig)

%ff stands for the difference of the two sides of the FOC of entrepreneurial utility
%maximization problem. Or the so-called efficiency condition of the
%entrepreneur
ff       =   (1-Gamma) * sp + ( sp*Gam_muG-1 )/( 1 - mu*omega*Fprime/(1-F) );
