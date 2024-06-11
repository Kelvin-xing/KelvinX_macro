function [diff,ks] = diseq_full(r,par,z,shocks)
global M_ oo_ 
% M_ and oo_ need to be in a global statement, because without it they are 
% not part of the memory of this function when Dynare is run below

k = zeros(par.T,1);

save parametervalues r par
%Saves parameter values to a file that will be read by Dynare
dynare Aiyagari_full noclearall
%Runs Dynare


%setting parameters for get_policy_rule_coefs.m
pert_order = 2;
VarsToUse = {...
     'k';...
     'c';...
      };
iprint = 0;

% The code in "get_policy_rule_coefs.m" will generate policy rule coefficients
% EXACTLY as Dynare writes them on the screen
% (unfortunately, the coefs in oo_ are a bit different than these)
% the matrix "decision" will contain the variables in the order specified here (so this
% may be different from the way it is shown by Dynare on the screen; up to you)

get_policy_rule_coefs
%Loads the decision rules to a matrix called 'decision'. 


%!!! IMPORTANT the starting point for the simulation must be something that
%is sensible for the 2nd-order perturbation solution you just calculated. 
%One sensible value is the steady state of the model you have specified in
%your Dynare file.
%The steady state specified in the mother program (which is for a different
%interest rate) could be very different and could very well be in the
%decreasing part of the second-order approximation. Then you'll get silly
%things.    
kss = decision(1)-decision(2);
k(1) = kss;
for t = 2:par.T
    k(t) = decision( 1,1) ...
         + decision( 3,1)*(k(t-1)-kss) ...
         + decision( 4,1)*(z(t-1)-par.z_ss)   ...
         + decision( 5,1)* shocks(t)   ...
         + decision( 6,1)*(k(t-1)-kss)^2 ...
         + decision( 7,1)*(k(t-1)-kss)*(z(t-1)-par.z_ss) ...
         + decision( 8,1)*(z(t-1)-par.z_ss)^2 ...
         + decision( 9,1)*(shocks(t))^2  ...
         + decision(10,1)*(k(t-1)-kss)*shocks(t) ...
         + decision(11,1)*(z(t-1)-par.z_ss)*shocks(t);
end
    
ks        = mean(k(par.T0:end));            % average capital supply (from HH problem)
kd        = (r/par.alpha)^(1/(par.alpha-1));    % implied capital
diff = ks-kd;

disp(diff)


