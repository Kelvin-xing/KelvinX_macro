%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
% This program looks at three models for equity prices when agents are
% boundedly rational
%
% Environment:
% - agents are risk neutral
% - dividends have (1) a persistent component that varies with the regime the
%   economy is in (low & high) and (2) standard AR(1) component with Normal
%   shocks
%
% Model I: Bayesian learning. 
% Agents know the correct formulas to determine equity prices, and the correct
% time-series model that generates the regimes, but they do not know the
% state the economy is in in a given period. After a regime change it will
% take some time before the agents have become more confident there has
% been a regime change. We will investigate how quickly this will happen
%
% Model II: "old-fashioned" least-squares learning. 
% Agents again know the
% correct pricing formula. Here, they are not aware of the dgp for dividends.
% In fact they erroneously believe that dividends will be constant and use a
% finite set of past observations to determine this "constant" level.
% The key aspect is that agents are ignorant of the regime changes. It
% would be easy to let them estimate a somewhat more complex dgp (say a
% standard AR(1) process), but this would not affect the point to be made.
%
% Model III: Adam-Marcet-Nicolini learning ("new" least-squares learning).
% Agents know the first-order condition, but not the correct formula for
% next-period's price. They use standard regression analysis. We assume
% that agents do know in which regime they are.
%
%
%% starting and setting parameter values
clc
clear all
close all

%parameter values
par.beta    = 0.99;                     % discount factor
par.rho     = 0.99;                     % probability of staying in regime
par.rho_z   = 0.5;                      % persistence of other stochastic component 
par.sigma   = 0.01;                     % volatility shock
par.mu_lo   =-0.01;                     % low state dividend parameter
par.mu_hi   =-par.mu_lo;                % high state dividend parameter

par.mu      = [par.mu_lo;par.mu_hi];

par.Tsimul  =   100;                     % time in regime is fixed in simulation (to understand what is going on)
                                         % if set equal to 1/(1-par.rho) then this is expected time the
                                         % stays in a regime (other values can be chosen)
par.T1      =   200;                     % par.T1 is one full cycle of a low and a high regime
                                         % AMN: fixed gain of 1/par.T1 
                                         % LS:  number of observations is par.T1/2. !!! Thus with LS
                                         % agents do not use a full cycle
par.T       = 10000;                     % total # of obs.

P           = [  par.rho 1-par.rho; ...
               1-par.rho   par.rho];     % transition matrix

%allocating memory          
dividend     = zeros(par.T,1);           % total dividends
mu_true      = zeros(par.T,1);           % value of mu (regime dependent dividend component)
z            = zeros(par.T,1);           % AR(1) component of dividends
price_ree    = zeros(par.T,1);           % REE price
price_lsq    = zeros(par.T,1);           % old-fashioned LS learning price
price_bay    = zeros(par.T,1);           % Bayesian learning price
price_amn    = zeros(par.T,1);           % AMN learning price
xi_tt        = zeros(par.T,2);           % beliefs about current state
xi_tnextt    = zeros(par.T,2);           % implied beliefs about next state
ind_mu       = zeros(par.T,1);           % -1,1 indicator for regime
IND_mu       = zeros(par.T,1);           %  1,2 indicator for regime

          
          
%% Solving for the REE
%
% the REE solution for p(t) is of the form B0_ree+B1_ree*z(t) 
% where B0_ree & B1_ree are 2x1 vectors. The top row correspond to the 1st
% regime and the bottom row to the second regime. (In other words there are
% two linear expressions, one for each regime)
% 

% XYZ Solve for B0_ree & B1_ree


B0_ree = ????;

B1_ree  = ????;

%% generate the random shocks

rng(20110810) % this is the new Matlab command to set the seed
% if you have an older version of Matlab you should use 
% randn('state',2011)

innovations  = randn(par.T,1)*par.sigma;

for t = 2:par.T
    z(t) = par.rho_z*z(t-1)+innovations(t);
end

% initializing variables for regime change
% the regime change is random but to better understand what is going on
% (and to get nice pictures) the time the system stays in each regime is
% fixed and set equal to its expected value
%
% !!! agents in the economy do not know we give them this specific set of
% realizations !!!

%% get ready for the simulation

itel       =  0;
ind_mu(1)  = -1;                                % ind_mu is -1 or 1
IND_mu(1)  = 1.5+0.5*ind_mu(1);                 % IND_mu is  0 or 1
                                                % 1 for  low state
                                                % 2 for high state

% the {-1,1} is convenient to generate the states, but the {1,2} is needed
% to indicate the corresponding element in the vector

% set initial beliefs about current (and implied next period) state:                                                
xi_tt(1,:)     = [0.5 0.5];
xi_tnextt(1,:) = [0.5 0.5];


%% the big simulation starts

for t = 2:par.T
    %% generate exogenous variables
    
    itel = itel + 1;
    if itel == par.Tsimul                       % after par.Tsimul periods itel is reset to 0
        itel = 0;
        ind_mu(t) = -ind_mu(t-1);
    else
        ind_mu(t) =  ind_mu(t-1);
    end
    IND_mu(t) = 1.5+0.5*ind_mu(t);              
    
    % the above makes sure that IND_mu(t) is 1 for par.Tsimul periods and
    % then 2 for par.Tsimul periods

    mu_true(t)     = par.mu(IND_mu(t));
    dividend(t)    = z(t) + mu_true(t);
    
    %% calculate rational expectations prices
    price_ree(t) = B0_ree(IND_mu(t)) + B1_ree(IND_mu(t))*z(t);
    price_amn(t)  = price_ree(t); % used for initialization only

    %% Bayesian learning
    
    % update probabilities for regime after observing 
    % current dividends and lagged z
    
    % XYZ: complete step I, II, and III
    % step I: calculate shocks conditional on possible regime
    
    res_lo         = ???;
    res_hi         = ???;
    
    % step II: update beliefs on current state and 
    % implied probs for next period's state
    
    xi_tt(t,:)     = ???;
    xi_tnextt(t,:) = P*xi_tt(t,:)';        
    
    % step III: calculate equilibrium price implied by these estimated probabilities
    % this is the actual price
    % we use the REE formulas for the two states weighted with perceived
    % probabilities of the state we are in
    
    price_bay(t)   = ??? ;

    %% with enough observations we can start doing the two types of least-squares learning
               
%     if t == par.T1
%         % observations up to par.T1 used to get initial estimate for
%         % precision matrix for AMN learning
%         X       = [ones(par.T1-1,1) ind_mu(1:par.T1-1) z(1:par.T1-1) z(1:par.T1-1).*ind_mu(1:par.T1-1)];        
%         RR      = X'*X;
%         % set initial beliefs about price process for AMN learning
%         beta_p  = [0; 0; 0; 0];
%     end
%     
%     if t > par.T1       
% 
%             %%  AMN learning
% 
%             % update pricing regression coefficients
%             
%             X        = [1 ind_mu(t-1) z(t-1) ind_mu(t-1)*z(t-1)];
%             Y        = price_amn(t-1);
%             note that beliefs in period t are based only on past
%             observations
%
%             XYZ complete the recursive least-squares formulas using fixed
%             % gain 1/par.T1
%             beta_p   = beta_p + ???;
%             RR       = RR     + ???;
%             
%             % calculate current price if next-period's price is believed to be determined by this process
%         
%             % no change in regime            
%             ind_mu_next = ind_mu(t);
%             IND_mu_next = 1.5+0.5*ind_mu_next;
%
%             XYZ: calculate the expected value of the RHS if there is no change in regime
%
%             RHS_next_1 = ???;
%
%
%             %    change in regime            
%             ind_mu_next = -ind_mu(t);
%             IND_mu_next = 1.5+0.5*ind_mu_next;
%
%             XYZ: calculate the expected value of the RHS if there is no change in regime
%
%             RHS_next_2 = ???;
%
%             % combine the two using the correct regime switching probabilities            
%             price_amn(t) = par.rho*RHS_next_1 + (1-par.rho)*RHS_next_2;
% 
%             %% OLD-Fashioned LS learning
%             div_mean    = mean(dividend(t-par.T1/2+1:t));
%             price_lsq(t) = par.beta*div_mean/(1-par.beta);
%             % !!! note that agents do not use enough obervations to cover
%             % both regimes (much more action this way)
%     end
end

%% make some figures

% figure(1)
% plot([price_ree(par.T1+1:end) price_amn(par.T1+1:end)])
figure(2)
plot([price_ree(par.T1+1:end) price_bay(par.T1+1:end)])
% figure(3)
% plot([price_ree(par.T1+1:end) price_lsq(par.T1+1:end)])

        
