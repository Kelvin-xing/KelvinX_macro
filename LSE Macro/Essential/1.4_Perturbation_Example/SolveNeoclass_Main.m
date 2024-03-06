%*************************************************************************************************************
%                                   `        Tools for Macroeconomists
%                                              PartI: The essentials
%                                               
%*************************************************************************************************************


% Parametrize and solve neoclassical growth model 

clear all
clc

%% 1. Parametrization
%=============================================================================================================

par.beta    = 0.99;         % discount factor
par.alpha   = 0.36;         % returns to scale in production
par.delta   = 0.025;        % depreciation rate
par.rhoz    = 0.95;         % autocorrelation of productivity shock
par.sigz    = 0.1;         % standard deviation of productivity shock
par.nu      = 1;            % relative risk aversion coefficient (1=log utility)

%% 2. Model solution and simulation
%=============================================================================================================

% steady state values
par.k       = ((1-par.beta*(1-par.delta)) ...   
                /(par.beta*par.alpha))^(1/(par.alpha-1));  % capital
par.c       = par.k^par.alpha - par.delta*par.k;                                    % consumption
par.z       = 1;

% return
% model solution
save params par
dynare neoclassModel.mod noclearall

return
oo_levels = oo_;

% model solution (variables in logs)
save params par
dynare neoclassModel_Logs.mod noclearall

oo_logs = oo_;

ytrue = (par.z+oo_levels.irfs.z_e).*(par.k+oo_levels.irfs.k_e).^par.alpha;

figure, plot(1:20,oo_levels.irfs.y_e/(par.k^par.alpha)*100,'k', ...
    1:20,oo_logs.irfs.y_e*100,'--r',1:20,(ytrue-par.k^par.alpha)/(par.k^par.alpha)*100,'-.b')
legend('approximation in levels','approximation in logs','using true production function')
ylabel('percent'),xlabel('time')
