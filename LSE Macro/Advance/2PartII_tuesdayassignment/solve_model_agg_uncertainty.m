%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
% Solving a model with both aggregate and invividual uncertainty using the
% algorithm of Krusell and Smith (1998).
%
%**************************************************************************
%
% notes: 
% 1) The problem of the individuals is solved using the dynare program "model_agg_uncertainty.mod" (second-order).
% 2) Simulated variables are all in absolute deviations from their own steady state levels.
% 3) To run this program, "disp_dr.m" is needed.
%
%**************************************************************************

warning off all
close all

% settings
T=10000;            % number of time periods in simulation step
N=1000;             % number of agents in simulation step
crit=0.000001;      % convergence criterion    

% parameter values
alpha = 0.25;
nu    = 1;
delta = 0.025;
beta  = 0.99;
rho   = 0.95;
k_ss  = (beta*alpha/(1-beta*(1-delta)))^(1/(1-alpha));
zeta0 = 0.01;
zeta1 = 0.3;
zeta2 = zeta1*exp(-zeta0*k_ss);
sig_e1   = 0.1;             
sig_e2   = 0.007;

rng(20100807,'philox')  % set the seed
%rng('default') % this sets the seed for the random number generator (to the default one when Matlab starts)


% initial values for coefficients law of motion for aggregate capital 
b_z=0.95;                            
b_0=1.7;
b_K=0.9;
save parametervalues k_ss beta nu delta zeta0 zeta1 zeta2 alpha rho b_0 b_K b_z sig_e1 sig_e2

coef_old = [b_0 b_K b_z]';

% reserve memory
k_sim_norm=zeros(N,T);                   % reserve memory for simulated individual capital stocks
k_sim_norm(:,1)=0;                  % k_sim_norm is normalized capital (see below); all agents start at the same level (doesn't matter much what it is)
Ka_sim=zeros(T,1);                  % reserve memory for simulated aggregate capital stock
Ka_sim(1)=mean(k_sim_norm(:,1));
z_sim=zeros(T,1);                   % reserve memory for aggregate productivity shock
z_sim(1)=0;
error=100;                          % initial value for error value in regression step

% draw shocks and simulate productivity
e1_sim=sig_e1*randn(N,T);               % draw innovations idiosyncratic productivity shock
e2_sim=sig_e2*randn(T,1);               % draw innovations aggregate productivity shock
for t=2:T
    % XXX
end

while error > crit
    % step 1: solve problem individual given a perceived law of motion for aggregate capital stock Ka
    save aggregatelaw b_0 b_K b_z sig_e1
    dynare model_agg_uncertainty.mod noclearall    % this generates the policy rule

    %setting parameters for get_policy_rule_coefs.m
    pert_order = 2;
    VarsToUse = {...
     'k';...
     'Ka';...
      };
    iprint = 0;

    % The code in "get_policy_rule_coefs.m" will generate policy rule coefficients
    % EXACTLY as Dynare writes them on the screen
    % (unfortunately, the coefs in oo_ are a bit different than these)
    % the matrix "decision" will contain the variables in the order specified here (so this
    % may be different from the way it is shown by Dynare on the screen; up to you)

    get_policy_rule_coefs
    %Loads the decision rules to a matrix called 'decision'. 
    %Alternatively, you use
    %
    % decision = get_policy_rule_coefs_fcn(pert_order,VarsToUse,iprint,M_,oo_)
    %
    % Here Dynare is called in the main program and decision is defined in
    % the main program so M_ and oo_ can be used to define decision
    % !! If Dynare is called _inside_ a function, then Matlab does not make
    % M_ and oo_ available inside function memory (only in main) so you
    % would have to add a global statement if that is what you want)

    
    
    % You can also get the coefficients of the decision rule yourself from what
    % Dynare has stored. But it is a bit tricky to figure out. For this
    % problem is would be as follows. (k is the second variable in the
    % DR-order which is something you cannot choose yourself so that is the
    % first thing you would have to check). But then for the constant the
    % ordering specified in the variable list is used.
    % k_i1 = 1;
    % k_i2 = 2;
    % DECISIONALT = [oo_.dr.ys(k_i1,1);oo_.dr.ghs2(k_i2,1)/2; oo_.dr.ghx(k_i2,:)'; oo_.dr.ghu(k_i2,:)'; oo_.dr.ghxx(k_i2,1)/2; oo_.dr.ghxx(k_i2,2); ... 
    %    oo_.dr.ghxx(k_i2,5)/2;oo_.dr.ghxx(k_i2,3);oo_.dr.ghxx(k_i2,6);oo_.dr.ghxx(k_i2,9)/2;oo_.dr.ghuu(k_i2,1)/2;oo_.dr.ghuu(k_i2,2);oo_.dr.ghuu(k_i2,4)/2;...
    %    oo_.dr.ghxu(k_i2,:)'];
    
    % step 2: simulate economy, consisting of N agents, over T periods
    for t=2:T
        % compute individual capital using decision rules from dynare
        %k_sim_norm is the normalized capital stock, that is, the actual
        %capital stock minus the deterministic steady state
        % the following line picks the right element according to the function defined above
    
    k_sim_norm(:,t)=[ones(N,1) k_sim_norm(:,t-1) ones(N,1)*z_sim(t-1) %XXX ]*decision(2:end,1);    
    
    % compute aggregate capital, taking account of the way dynare writes
    % decision rules.
    Ka_sim(t)= %XXX
    end    
    
    % step 3: run a regression and update law of motion for Ka
    % XXX
    
    % step 4: evaluate convergence
    
    error = %XXX
    disp(error)
    disp('  ')
    disp('program has paused')
    disp('hit any key to continue, but check whether error is getting smaller')
    disp('  ')
    pause

end

             