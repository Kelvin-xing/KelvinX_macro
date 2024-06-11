%*************************************************************************************************************
%                                   `        Tools for Macroeconomists
%                                              PartI: The essentials
%                                               
%*************************************************************************************************************


% Parametrize and solve neoclassical growth model - simulate data and use simulated data to estimate some of 
% the parameters - using ML or Bayesian estimation

clear all
clc

%% 1. Parametrization
%=============================================================================================================

par.beta    = 0.99;         % discount factor
par.alpha   = 0.36;         % returns to scale in production
par.delta   = 0.025;        % depreciation rate
par.rhoz    = 0.95;         % autocorrelation of productivity shock
par.sigz    = 0.01;         % standard deviation of productivity shock
par.nu      = 1;            % relative risk aversion coefficient (1=log utility)

%% 2. Model solution and simulation
%=============================================================================================================

% steady state values
par.k       = ((1-par.beta*(1-par.delta)) ...   
                /(par.beta*par.alpha))^(1/(par.alpha-1));  % capital
par.c       = par.k^par.alpha - par.delta*par.k;                                    % consumption
par.z       = 1;

% model solution
save params par
dynare neoclassModel.mod noclearall

% simulated time-series
c       = oo_.endo_simul(1,1000:1265);
k       = oo_.endo_simul(2,1000:1265);
y       = oo_.endo_simul(3,1000:1265);

pause
%% 3. Estimation using the output time-series via Maximum Likelihood
%=============================================================================================================

% 3.1 Estimating only sigma_z
%------------------------------
 
% save data file for loading within dynare
save y y 

save params par
dynare neoclassModel_estim.mod noclearall

% return

% 3.2 Estimating sigma_z, rho_z, delta, alpha
%----------------------------------------------
pause

% save data file for loading within dynare
save y y 

save params par
dynare neoclassModel_estim2.mod noclearall

