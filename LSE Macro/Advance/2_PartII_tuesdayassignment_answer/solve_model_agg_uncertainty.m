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
clear
warning off all
close all
 
% settings
T        = 10000; % number of time periods in simulation step
N        = 1000;  % number of agents in simulation step
crit     = 1e-6;  % convergence criterion    
error    = 100;   % initial value for error value in regression step
Tdiscard = 100;    % discard first X periods for the regeression
adj      = 0.25;  % smoothing parameter 

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

z_ss=1;

% initial values for coefficients law of motion for aggregate capital 
b_z=0.95;                            
b_0=1.7;
b_K=0.9;

coef_history=[b_0 b_K b_z];

save parametervalues k_ss beta nu delta zeta0 zeta1 zeta2 alpha rho b_0 b_K b_z sig_e1 sig_e2
rng(20100807,'philox')  % set the seed
%rng('default') % this sets the seed for the random number generator (to the default one when Matlab starts)

coef_old = [b_0 b_K b_z]';
new_coef=coef_old;

% reserve memory
k_sim_norm  = zeros(N,T);               % reserve memory for simulated individual capital stocks
Ka_sim_norm = mean(k_sim_norm);         % innitialise aggregate capital
z_sim       = zeros(1,T);               % reserve memory for aggregate productivity shock
z_sim(1)    = 0;

% draw shocks and simulate productivity
e1_sim=sig_e1*randn(N,T);               % draw innovations idiosyncratic productivity shock
e2_sim=sig_e2*randn(1,T);               % draw innovations aggregate productivity shock
for t=2:T
    z_sim(t) = z_ss*(1-rho) + rho*z_sim(t-1)+e2_sim(t);
end

z_sim_norm=z_sim-z_ss;

%state vector at time t consiste of both aggregate variables (KA,e1,e2,z) and idiosyncratic k

% the following line defines the function. At the end of this program
% you'll find a simple example explaining this structure

state_vector = @(t,k,z_sim,Ka,e1,e2) ...
                    [ones(N,1), ...
                        k(:,t-1), ...
                        ones(N,1)*z_sim(t-1), ... 
                        ones(N,1)*Ka(t-1), ...
                        e1(:,t), ... 
                        ones(N,1)*e2(t), ...
                        k(:,t-1).^2, ...
                        k(:,t-1)*z_sim(t-1), ...
                        ones(N,1)*z_sim(t-1)^2, ...
                        k(:,t-1)*Ka(t-1), ... 
                        ones(N,1)*z_sim(t-1)*Ka(t-1), ...
                        ones(N,1)*Ka(t-1)^2, ...
                        e1(:,t).^2, ...
                        e1(:,t)*e2(t), ...
                        ones(N,1)*e2(t)^2, ...
                        k(:,t-1).*e1(:,t), ... 
                        k(:,t-1)*e2(t), ...
                        z_sim(t-1)*e1(:,t), ...
                        ones(N,1)*z_sim(t-1)*e2(t), ...
                        Ka(t-1)*e1(:,t), ... 
                        ones(N,1)*Ka(t-1)*e2(t)...
                     ];

while error > crit        
    % step 1: solve problem individual given a perceived law of motion for aggregate capital stock Ka
    save aggregatelaw b_0 b_K b_z
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
        k_sim_norm(:,t) = state_vector(t,k_sim_norm,z_sim_norm,Ka_sim_norm,e1_sim,e2_sim)*decision(2:end,1);       
        Ka_sim_norm(t) = mean(k_sim_norm(:,t)); %update the value of the per capita capital stock
%         disp(Ka_sim_norm(t))
%         pause
    end
  
    % compute aggregate capital, taking account of the way dynare writes decision rules.        
    k_ss  = decision(1,1)-decision(2,1); 
    k_sim = k_sim_norm+k_ss;
    
    Ka_sim = mean(k_sim,1);    
    Ka_ss = b_0/(1-b_K);
 %   Ka_sim_norm = mean(k_sim_norm,1);

    % step 3: run a regression and update law of motion for Ka
    X = [ones(1,T-1-Tdiscard)', Ka_sim(1+Tdiscard:end-1)', z_sim(2+Tdiscard:end)'-1];
    Y = Ka_sim(2+Tdiscard:end)';
    
    new_coef = regress(Y,X);
    b_0      = adj*new_coef(1)+(1-adj)*b_0; 
    b_K      = adj*new_coef(2)+(1-adj)*b_K;
    b_z      = adj*new_coef(3)+(1-adj)*b_z;
    
    coef_history = [coef_history;b_0,b_K,b_z];
    
    % evaluate convergence
    error = norm(coef_history(end-1,:)-[b_0,b_K,b_z]);
    disp(error)
    disp('  ')
    disp('program has paused')
    disp('hit any key to continue, but check whether error is getting smaller')
    disp('  ')
    pause
end

disp([b_0 b_K b_z])
dynare model_agg_uncertainty_IRF.mod noclearall


%% simple example explaining the function structure used above

check = zeros(5,3,8);

N=5;
k=randn(5,7);
z=randn(1,7);

% the folllowing defines the function state_vector which says that
% state_vector(t,k,z) is a matrix with
% an Nx1 column of ones in the first column
% an Nx1 column with the double lagged value of k(:,t)
% an Nx1 column with the leaded value of z(1,t)

state_vector = @(t,k,z) ...
                    [ones(N,1),  ...
                        k(:,t-2), ...
                        ones(N,1)*z(:,t+1)];
                 
for t = 3:6
    check(:,:,t) = state_vector(t,k,z);
end
disp(k)
disp(z)
 disp(check(:,:,5))


% 
             