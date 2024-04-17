function [ys,check]=bgg_rbc_steadystate(ys,exe)

% compute s.s. This is the steady state file for the mod file;
% this file will automatically invoked by dynare and then s.s. 
% are returned.

global M_ options_
check = 0;

%the following statements ensure that the variables are not interpreted by
%MATLAB as functions... they will be reassigned values
%you do not need worry about values here.
beta=5;
alpha=4;
gamma=4; 
mu=1;
% Here we load the values of the deep parameters in a loop
% from dynare pre-processing info M_.
Np = M_.param_nbr;                                            
for i = 1:Np
    paramname = deblank(M_.param_names(i,:));
    eval([ paramname ' = M_.params(' int2str(i) ');']);
end

z=1;
R=(1/beta)-1;
[Rk,omegabar,G,F,Gamma,Gam_muG,Fprime,kbar,n,c] = steadystate(1+R,sigma_ss,mu,alpha,gamma,delta,z);
Rk=Rk-1;
i=delta*kbar;
sigma=sigma_ss;
credit = kbar-n;
spread = ((1+Rk)/(1+R))*(kbar/credit)*omegabar;
bankrupt=F;
GDP=c+i;
wedge=(1-(1+R)/(1+Rk));


% Define the steady state values of the endogenous variables of the model.
%@2013-9-2 This steady state file will automatically invoked by Dynare, ys
%are returned to the Dynare.
Ne = M_.orig_endo_nbr;
ys = zeros(Ne,1);
% endoleadcount = 2;
% nonauxcount = 0;
for indexvar = 1:Ne
    varname = deblank(M_.endo_names(indexvar,:));
        eval(['ys(' int2str(indexvar) ') = ' varname ';']);
end