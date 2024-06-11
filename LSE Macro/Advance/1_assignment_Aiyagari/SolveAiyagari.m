%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%==========================================================================
%                  Solving an Aiyagari type model
%==========================================================================

% This program solves an Aiyagari type model (idiosyncratic shocks and a 
% penalty function on capital accumulation). 
%
% First an increase in volatility is considered in the equilbrium model.
%
% Second, using the equilibrium interest rate of the low-volatility
% economy, we then resolve the two economies. Thus, we now measure the
% partial equilibrium response of the increase in volatility

clear
clc
% M_ and oo_ need to be added to global memory. Without doing this, their values would not be known
% inside the function diseq_full.m even though Dynare is run inside that function
global M_ oo_

%--------------------------------------------------------------------------
% 1. Parameters
%--------------------------------------------------------------------------

par.T  = 100000;                        % length of simulation
par.T0 = 10001;                         % start of sample

par.tol   = 1E-4;                       % tolerance level for convergence

par.beta  = 0.99;                       % discount factor
par.alpha = 0.25;                       % curvature prod. function
par.delta = 0.025;                      % depreciation rate
par.rho   = 0;                          % persistence parameter on shocks
par.nu    = 0.5;                        % coefficient of relative risk aversion

k_ss      = (par.beta*par.alpha/(1-par.beta*(1-par.delta)))^(1/(1-par.alpha));  %steady state capital 
par.kini  = k_ss;                       %initial value used in *.mod file

par.z_ss  = 1;                          % steady state value of shocks

par.zeta0 = 0.10;                       % parameter penalty function
par.zeta1 = 0.003;                      % parameter penalty function
par.zeta2 = par.zeta1*exp(-par.zeta0*k_ss);     % parameter penalty function

% zeta2 is chosen such that steady state value of capital is not affected

par.r_ss = (1-par.beta*(1-par.delta))/par.beta;  % standard steady state value
r        = 0.0351;                               % initial value for r

% allocating memory and shock values
z     = ones(par.T,1);
k     = zeros(par.T,1);

ks_partial = zeros(2,1);
ks_ge      = zeros(2,1);
r_ge       = zeros(2,1);

sigs = [0.001;0.3];             % values for standard deviations of shocks

%--------------------------------------------------------------------------
% 2. Generate innovations
%--------------------------------------------------------------------------

rng(20100807,'philox') % set the seed of the random number generator
innovations = randn(par.T,1);     % defining innovation series

%--------------------------------------------------------------------------
% 3. Loop over sigma for equilibrium r (GENERAL EQUILIBRIUM)
%--------------------------------------------------------------------------
% This section solves for the general equilibrium interest rate. It
% contains an algorithm that searches for r, at which capital supply from
% the household is the same as capital demand of the representative firm.
% This is done for two values of sigma (the idiosyncratic risk).
%--------------------------------------------------------------------------

for i = 1:2
    
% defining shocks
par.sigshock = sigs(i);                 % choosing standard deviation
shocks   = par.sigshock*innovations;    
for t = 2:par.T
    z(t) = 1-par.rho + par.rho*z(t-1) + shocks(t);
end
%    it often is good to adjust z such sample mean is correct. Thus 
%    z = z + 1 - mean(z); 
%    but then you have to use only z to construct choices not z(-1) and
%    shocks as is done here

%--------------------------------------------------------------------------
% Solve for equilibrium r
%--------------------------------------------------------------------------

r_l       = r-0.001;    % lower bound for r at which ks < kd
r_u       = r+0.001;    % upper bound for r at which ks > kd
err       = 100;        % initial error
r_x       = r_u;

[diseq_r_l,ks_l] = diseq(r_l,par,z,shocks);
[diseq_r_u,ks_u] = diseq(r_u,par,z,shocks);
[diseq_r_x,ks_x] = diseq(r_x,par,z,shocks);
prodct    = diseq_r_l*diseq_r_u;


if prodct > 0 % A warning message if initial interval for r is too narrow
    disp(' root not in between r_l and r_u')
    disp(' adjust r_l or r_u')
    stop
end

    disp(' ')
    disp([r_l r_x r_u diseq_r_l diseq_r_x diseq_r_u])
    %pause(2)

while err > par.tol
        
    r_x_old = r_x;
    r_x     = r_u - diseq_r_u*(r_u-r_l)/(diseq_r_u-diseq_r_l);
    [diseq_r_x,ks_x] = diseq(r_x,par,z,shocks);

    disp(' ')
    disp(diseq_r_x)
    disp([r_l r_x r_u diseq_r_l diseq_r_x diseq_r_u])
    %pause(2)
    
    if r_x ~= 0
        err = abs((r_x - r_x_old)/r_x) * 100;
    else
        err = abs(r_x - r_x_old);
    end
        
    test = diseq_r_l*diseq_r_x;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below fill in what r_x should be if test<0 and what it should be if
% test>0. Do the same for diseq_r_x
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if test == 0
        err = 0;  % root is at r_x
    elseif test < 0
         = r_x;        % root is below r_x
         = diseq_r_x;
    else
         = r_x;        % root is above r_x
         = diseq_r_x;
    end

end

ks_ge(i) = ks_x; % general equilibrium capital
r_ge(i)  = r_x;% general equilibrium interest rate
end

% save values for proj_PE.m This just saves the Dynare's results. They
% will be used later as starting values for projection.
save r_ge
save ks_ge
save par

%--------------------------------------------------------------------------
% 4. Loop over sigma for fixed r (PARTIAL EQUILIBRIUM)
%--------------------------------------------------------------------------
% Here we solve for the partial equilibrium twice. First, we take the
% aggregate interest rate from the model with low uncertainty and solve for
% the individual problem. Second, we increase the uncertainty, but keep the
% level of r fixed at the equilibrium value from the low uncertainty case
% (i.e., solve for the partial equilibrium).
%--------------------------------------------------------------------------

% r = r_ge(1);        % fixing interest rate
% 
% for i = 1:2
% 
% % i = 1: we use r=r_ge(1), i.e., the low volatility equilibrium interest
% %        rate, so we simply resolve at equilibrium value
% % i = 2: we increase uncertainty but do not adjust r, i.e. partial eq.
% 
%     par.sigshock = ;
%     shocks   = ;
%     
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % Below solve the model with dynare for the fixed interest rate and the given
% % standard deviation of the shocks
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
%     save
%     %Saves parameter values to a file that will be read by Dynare
%     dynare
%     %Runs Dynare

%setting parameters for get_policy_rule_coefs.m
%pert_order = 2;
%VarsToUse = {...
%     'k';...
%     'c';...
%      };
%iprint = 0;

% % The code in "get_policy_rule_coefs.m" will generate policy rule coefficients
% % EXACTLY as Dynare writes them on the screen
% % (unfortunately, the coefs in oo_ are a bit different than these)
% % the matrix "decision" will contain the variables in the order specified here (so this
% % may be different from the way it is shown by Dynare on the screen; up to you)

%get_policy_rule_coefs

% %Loads the decision rules to a matrix called 'decision'. 

% 
%     for t = 2:par.T
%     z(t) = 1 - par.rho + par.rho*z(t-1) + shocks(t);
%     end
%     z = z + 1 - mean(z); %impose mean is not subject to sampling variation
% 
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % Below specify the steady state capital stock and simulate the economy.
% % That is, write a recursion that will generate a long time series for
% % capital.
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     
%     kss = ;
%     k(1) = kss;
%     for t = 2:par.T
% 
%     end
% 
%     ks_partial(i) = mean(k(par.T0:end));  % average capital supply (from HH problem)
% 
% end
% 
% disp(' ')
% disp('Partial equilibrium capital supply')
% disp([ks_partial(1), 100*(ks_partial(2)-ks_partial(1))/ks_partial(1)])
% disp('General equilibrium capital supply')
% disp([ks_ge(1), 100*(ks_ge(2)-ks_ge(1))/ks_ge(1)])

