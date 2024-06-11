%%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
% 
% sunspotlearning_sinus_linear1.m
%
% this program uses adaptive learning to see whether a sunspot can be
% learned. The model is from "Indeterminacy and E-stability in Real
% Business Cycle Models with Factor-Generated Externalities by
% Bruce McGough, Qinglai Meng, and Jianpo Xue who show the linearized
% version has learnable sunspot solutions
%
%% overview program
%
%  1. set parameter values algorithm
%  2. set parameter values model & steady state
%  3. get the linearized system
%  4. generate time series with linearized system
%  5. check accuracy linearlized system (optional)
%  6. get the nonlinear PEA approximation
%  7. check accuracy of PEA approximation


%% 
%
%   1. Parameter values algorithm
%

clc
clear all
close all

initialvalues        = 'linear'; % use either linear or previous
version              = 'case1';
% case1:  first-order for non-sunspot &  first-order for sunspot part of approximation & MMX timing
%         discrete shocks
% case2:  first-order for non-sunspot &  first-order for sunspot part of approximation & MMX timing
%         shocks with continuous support
% case3:  first-order for non-sunspot &  first-order for sunspot part of approximation & standard  timing
%         shocks with continuous support
% case4:  higher-order for non-sunspot & higher-order for sunspot part of approximation & standard timing
%         shocks with continuous support
accuracylinear       = 'no';
accuracyprojection   = 'no';
% yes/no for Dynamic Euler equation test

T              = 11000;
T1             =  1000;
conv_criterion =  1E-4;
dampening      =  0.5;
par.N_herm     =  7; % number of quadrature nodes for accuracy test
                
options = optimset('Display','off','MaxFunEvals',1E5,'MaxIter',...
                  1E5,'TolFun',1E-6,'TolX',1E-6);              
              
% these are parameters for the minimization routine
% this is a numerical problem within the bigger numerical problem. 
% such innerloop problems should be done more accurately
%
%% 
%  2. Parameter values model
%

switch version
    case 'case1'
        par.eta_dzeta  = -0.00375;         % maximum impact sunspot on MU
    case 'case2'
        par.eta_dzeta  = -0.00375;         % maximum impact sunspot on MU
    case 'case3'
        par.eta_dzeta  = -0.00375;         % maximum impact sunspot on MU
    case 'case4'
        par.eta_dzeta  = -0.00375;         % maximum impact sunspot on MU
end

par.a          =  0.3;              % capital share production        
par.b          =  0.7;              % labor share production
par.delta      =  0.025;            % depreciation
par.dfactor    =  0.99;             % discount factor
par.beta       =  1.1;              % aggregate labor coefficient
par.chi        =  0.08;             % disutility of labor
        
% the following two parameters are adjusted in MMX
        
par.alpha      =  0.25;             % aggregate capital coefficient        
par.nu         =  0.1;              % risk aversion (sigma in MMX)

par.vola       =  1;             % standard deviation of the shock
par.lin_dzeta  =  0.02;          % impact sunspot on cons in linearized ss

% labda1 & labda2 are scalingscoefficients to get L_ss = K_ss = 1
% note that these scalingscoefficients only affect the steady state
% and not of the dynamic properties of the model
        
par.L_ss       = 1;
par.K_ss       = 1;
par.labda1     = (1-par.dfactor*(1-par.delta))/(par.dfactor*par.a);
par.Y_ss       = par.labda1;
par.C_ss       = par.Y_ss - par.delta*par.K_ss;
par.labda2     = par.b*par.labda1*par.C_ss^(-par.nu);
        
% check steady state (these numbers should be 0)
xxx = [ ...
 -par.C_ss^(-par.nu)+par.dfactor*par.C_ss^(-par.nu)*(par.a*par.labda1+1-par.delta),  ...
 -par.labda2+par.b*par.labda1*par.C_ss^(-par.nu)                                         ,  ...
  par.labda1+(1-par.delta)-par.C_ss-1                                         ];
xxx0 = sprintf('model eqs. at steady state: %10.8f %10.8f %10.8f',xxx);
disp(' ');
disp(xxx0)
disp(' ');
        
 %%
 %   3. Get log-linearized system
 %
 
% linearize budget constraint: k(t+1) = dd_c*c(t) + dd_k*k(t) + dd_l*l(t)
dd_c = -par.C_ss/par.K_ss;
dd_k = (1-par.delta)+par.alpha*par.Y_ss/par.K_ss;
dd_l =              +par.beta *par.Y_ss/par.K_ss;

% linearize FOC labor: l(t) = ff_k*k(t) + ff_c*c(t)
ff_k =  par.alpha/(par.chi+1-par.beta);
ff_c = -par.nu   /(par.chi+1-par.beta);

% substitute out labor from budget constraint:
d_c  = dd_c + dd_l*ff_c;
d_k  = dd_k + dd_l*ff_k;

%linearize Euler equation: c(t) = bb_c*c(t+1) + bb_k*k(t+1) + bb_l*l(t+1)
bb_c =  1;
bb_k = -(par.alpha-1)*par.dfactor*par.a*(par.Y_ss/par.K_ss)/par.nu;
bb_l =  -par.beta    *par.dfactor*par.a*(par.Y_ss/par.K_ss)/par.nu;

% substitute out labor from Euler equation:
b_c  = bb_c + bb_l*ff_c;
b_k  = bb_k + bb_l*ff_k;

% combine Euler equation and budget constraint in matrix equation:
LL_next = [   1, 0  ; b_k, b_c];
LL_now  = [ d_k, d_c;   0, 1  ];

% [k(t+1) c(t+1)]' = JJ*[k(t) c(t)]' + [0 sunspot(t+1)]'
JJ = LL_next\LL_now;

% What is the solution to the model?
% 
% 1. If JJ has two eigen values with modulus less than 1, then you are done
% because this would be the solution. 
% In this case there is indeterminacy. Not only can you add a 
% sunspot shock, you are also free to pick the initial value for
% consumption arbitrarily. This is the case we focus on in this program
%
% 2. If JJ has one eigen value with modulus bigger, then you get the more
% standard case of a unique solution but you would have to do a bit more
% work. Let
% JJ = P*LABDA*Q, where
% LABDA is the diagonal matrix with eigenvalues
% P is the matrix with associated eigenvectors
% Q is the inverse of P
%
% to avoid explosive solution you have to 
% (i) set sunspot(t)=0 for all t
% (ii) set Q1*[k(t) c(t)]'=0, where Q1 is the row of Q corresponding to the
% eigen value > 1. That is, you are no longer free to pick any arbitrary
% c(t) but have to pick c(t) as a function of k(t).
 
%%
%  4. Generate time series using linearized solution
%

 x_tilde = zeros(T,2);
 X       = zeros(T,2);

% generate sunspot variable

rng(20110108) % this replace the former Matlab command to set the seed
% if you have an older version of Matlab replace it with 
% randn('seed', 20110108);

%dzeta = -1+2*(rand(T,1)-0.5>0); % to generate RV with discrete support
%discrete support can only be used for first-order approximation

switch version
    case  'case1'
        dzeta = -1+2*(rand(T,1)-0.5>0); % to generate RV with discrete support
    otherwise
        dzeta = par.vola*randn(T,1);       % to generate RV with continuous support
end

% generate time paths

for t = 2:T
    x_tilde(t,:) = (JJ*x_tilde(t-1,:)' + [0;par.lin_dzeta*dzeta(t)])';
end

X = [x_tilde(:,1) + log(par.K_ss), ... 
     x_tilde(:,2) + log(par.C_ss)];
X =  exp(X);

%%
%  5. Dynamic Euler equation test for linearlized system
%

switch accuracylinear
    case 'yes'
    
        par.JJ   = JJ;
        par.ff_k = ff_k;
        par.ff_c = ff_c;

        C_simulB = zeros(T  ,1);
        L_simulB = zeros(T  ,1);
        K_simulB = zeros(T+1,1);

        K_simulB(1) = par.K_ss;
        K_simulB(2) = par.K_ss;
        C_simulB(1) = par.C_ss;
     
        for t = 2:T;
        %generate C(t), K(t), & L(t):
        
            k_tilde       = par.JJ(1,1)*log(K_simulB(t-1)/par.K_ss)+par.JJ(1,2)*log(C_simulB(t-1)/par.C_ss);
            c_tilde       = par.JJ(2,1)*log(K_simulB(t-1)/par.K_ss)+par.JJ(2,2)*log(C_simulB(t-1)/par.C_ss)+par.lin_dzeta*dzeta(t);
        
            K_t   = exp(k_tilde + log(par.K_ss));
            C_t   = exp(c_tilde + log(par.C_ss));

            par.C_t     = C_t;
            par.K_t     = K_t;
            par.K_next  = exp(log(par.K_ss)+par.JJ(1,1)*log(K_t/par.K_ss)+par.JJ(1,2)*log(C_t/par.C_ss));
            
            switch version
                case 1
                    CondExp     = numi_case1(@RHSrealization_sunspot_lin,par);            
                otherwise
                    CondExp     = numi(@RHSrealization_sunspot_lin,par);            
            end
        
            C_simulB(t)  = CondExp^(-1/par.nu);

            L_simulB(t)  = (par.b*CondExp*par.labda1*K_simulB(t)^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
            K_simulB(t+1)= par.labda1*K_simulB(t)^par.alpha*L_simulB(t)^par.beta+(1-par.delta)*K_simulB(t)-C_simulB(t);
%           disp([100*abs(log(par.K_next/K_next_alt)) 100*abs(log(C_simulB(t)/C_alt))])
            %disp([X(t,1) K_simulB(t) X(t,2) C_simulB(t)])
            %pause
        end
        figure(1)
        plot([X(1001:1500,1) K_simulB(1001:1500)])
        figure(2)
        plot([X(1001:1500,2) C_simulB(1001:1500)])
        mean_Klinear = 100*mean(abs(log(X(:,1)./K_simulB(1:end-1))));
        mean_Clinear = 100*mean(abs(log(X(:,2)./C_simulB(1:end  ))));
         max_Klinear = 100* max(abs(log(X(:,1)./K_simulB(1:end-1))));
         max_Clinear = 100* max(abs(log(X(:,2)./C_simulB(1:end  ))));
         disp(' ')
         xxx0 = sprintf(' standard deviation of K & C:    %7.3f %7.3f',std(log(X(:,1))),std(log(X(:,2))));
         xxx1 = sprintf('  mean perc. error for K & C:    %7.3f %7.3f',mean_Klinear,mean_Clinear);
         xxx2 = sprintf('   max perc. error for K & C:    %7.3f %7.3f', max_Klinear, max_Clinear);
         disp(xxx0)
         disp(xxx1)
         disp(xxx2)
         disp(' ')
         disp('paused - hit any key to continue')
        pause
end

disp(' ')


%%
%   6. PEA
%

% initial work for PEA

switch initialvalues
    case 'linear' 
        % use the linearized solution (for log(C(t)) to get initial conditions for the policy
        % function of period t marginal utility:
        
        %common coefficients for non-sun-spot part:
        eta1 = [-par.nu*log(par.C_ss);-par.nu*JJ(2,1);-par.nu*JJ(2,2)];
        switch version
            case 'case1'
                %coefficients for sun-spot part:
                eta2 = 0.9;
            case 'case2'
                %coefficients for sun-spot part:
                eta2 = 0.9;
            case 'case3'
                %coefficients for sun-spot part:
                eta2 = 1.8;
            case 'case4'
                eta1 = [eta1;0;0;0];
                %coefficients for sun-spot part:
                eta2 = [0.1; 0; 0];
        end
        %all coefficients5
        eta = [eta1;eta2];        
    case 'previous'
        switch version
            case 'case1'
                load eta_sinus1
            case 'case2'
                load eta_sinus2
            case 'case3'
                load eta_sinus3
            case 'case4'
                load eta_sinus4
        end
end

C_simul = zeros(T  ,1);
L_simul = zeros(T  ,1);
Y_simul = zeros(T  ,1);
K_simul = zeros(T+1,1);
K_temp  = zeros(T+1,1);

K_simul(1) = par.K_ss;
K_simul(2) = par.K_ss;
C_simul(1) = par.C_ss;

 error_iter = 1000;

%
% start the big iteration (hopefully till convergence)
%

while error_iter > conv_criterion

    % step I of each iteration: generate time series
    
    for t = 2:T;
        
        % define state used to make forecast:
        % (note that MMX use k(t-1) not k(t))
        
        switch version
            case 'case1'
                state1        = [ 1 log(K_simul(t-1)/par.K_ss) log(C_simul(t-1)/par.C_ss)     ];
                state2        =   dzeta(t)*pi/2 ;
            case 'case2'
                state1        = [ 1 log(K_simul(t-1)/par.K_ss) log(C_simul(t-1)/par.C_ss)     ];
                state2        =   dzeta(t)*pi/2 ;
            case 'case3'
                state1        = [ 1 log(K_simul(t)  /par.K_ss) log(C_simul(t-1)/par.C_ss)     ];
                state2        =   dzeta(t)*pi/2 ;
            case 'case4'
                state1        = [ 1 log(K_simul(t)  /par.K_ss) log(C_simul(t-1)/par.C_ss) ...
                                   (log(K_simul(t)  /par.K_ss))^2                         ...
                                   (log(C_simul(t-1)/par.C_ss))^2                         ...
                                    log(C_simul(t-1)/par.C_ss)*log(K_simul(t)  /par.K_ss)     ];
                state2        = [ dzeta(t)*pi/2 dzeta(t)^2 dzeta(t)*log(K_simul(t)/par.K_ss)];
        end
        
        N_s1          = size(state1,2);
        N_s2          = size(state2,2);
        N_s           = N_s1+N_s2;
        mu_t          = exp(state1*eta1)+par.eta_dzeta*sin(state2*eta2);
        C_simul(t)   = (mu_t)^(-1/par.nu);
%        the following four lines give solutions for K & L according to linearized solution         
%        k_tilde      = d_k *log(K_simul(t-1)/par.K_ss)+d_c *log(C_simul(t-1)/par.C_ss);
%        K_temp(t)    = exp(log(par.K_ss) + k_tilde);
%        l_tilde    = ff_k*log(K_simul(t)  /par.K_ss)+ff_c*log(C_simul(t)  /par.C_ss);
%        L_simul(t) = exp(log(par.L_ss) + l_tilde);

        L_simul(t) = (par.b*C_simul(t)^(-par.nu)*par.labda1*K_simul(t)^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
        K_simul(t+1) = par.labda1*K_simul(t)^par.alpha*L_simul(t)^par.beta+(1-par.delta)*K_simul(t)-C_simul(t);
        
        % finally generate the thing inside the conditional expectation:
        Y_simul(t) = par.dfactor*C_simul(t)^(-par.nu)*(par.a*par.labda1*K_simul(t)^(par.alpha-1)*L_simul(t)^par.beta+1-par.delta);
        %disp([C_simul(t) X(t,2) dzeta(t)])
    end
    
    % step II of each iteration: update beliefs by running NLLS regression
    
    % give initial values for 3 free values of eta and set fixed parameter
    eta_in     = eta;

    Y_simul  = Y_simul(T1:T); 
    
    switch version
        case 'case1'
            X_simul1 = [ones(T-T1+1,1) log(K_simul(T1-2:T-2)/par.K_ss) log(C_simul(T1-2:T-2)/par.C_ss)          ];
            X_simul2 =  dzeta(T1-1:T-1)*pi/2;
        case 'case2'
            X_simul1 = [ones(T-T1+1,1) log(K_simul(T1-2:T-2)/par.K_ss) log(C_simul(T1-2:T-2)/par.C_ss)          ];
            X_simul2 =  dzeta(T1-1:T-1)*pi/2;
        case 'case3'
            X_simul1 = [ones(T-T1+1,1) log(K_simul(T1-1:T-1)/par.K_ss) log(C_simul(T1-2:T-2)/par.C_ss)          ];
            X_simul2 =  dzeta(T1-1:T-1)*pi/2;
        case 'case4'
            X_simul1 = [ones(T-T1+1,1) log(K_simul(T1-1:T-1)/par.K_ss) log(C_simul(T1-2:T-2)/par.C_ss) ... 
                                      (log(K_simul(T1-1:T-1)/par.K_ss)).^2                             ...
                                      (log(C_simul(T1-2:T-2)/par.C_ss)).^2                             ...
                                       log(C_simul(T1-2:T-2)/par.C_ss).*log(K_simul(T1-1:T-1)/par.K_ss)         ];
            X_simul2 = [dzeta(T1-1:T-1)*pi/2 dzeta(T1-1:T-1).^2 dzeta(T1-1:T-1).*log(K_simul(T1-1:T-1)/par.K_ss)];
    end
    
    % inputs of rssfcn are eta,Y_simul,X_simul,X_dzeta, and par.eta_dzeta and
    % minimization is with respect to eta
    eta_hat    = fminsearch(@(eta) rssfcn_sinus(eta,par.eta_dzeta,Y_simul,X_simul1,X_simul2,N_s1,N_s2),eta_in,options);    

    error_iter = sum(abs(eta_hat-eta));
    xxx1 = sprintf('error     = %4e',error_iter);
    xxx2 = sprintf('eta       = %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f',eta);  
    xxx3 = sprintf('eta_hat   = %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f',eta_hat);  
    disp(xxx2)
    disp(xxx3)
    disp(xxx1)
    %pause
    eta   = (1-dampening)*eta_hat + dampening*eta;
    eta1 = eta(     1:N_s1     );
    eta2 = eta(N_s1+1:N_s1+N_s2);

 end
 
switch version
    case 'case1'
        save eta_sinus1 eta eta1 eta2
    case 'case2'
        save eta_sinus2 eta eta1 eta2
    case 'case3'
        save eta_sinus3 eta eta1 eta2
    case 'case4'
        save eta_sinus5 eta eta1 eta2
end

%% 
%   7. Do the dynamic Euler equation test
%

switch accuracyprojection
    case 'yes'

        par.eta1 = eta1;
        par.eta2 = eta2;

        %dzeta =  par.vola*randn(T,1); % get fresh draw for more powerful test

        C_simulA = zeros(T  ,1);
        L_simulA = zeros(T  ,1);
        K_simulA = zeros(T+1,1);

        K_simulA(1) = par.K_ss;
        K_simulA(2) = par.K_ss;
        C_simulA(1) = par.C_ss;
     
        for t = 2:T            
            %generate C(t), K(t), & L(t):
            
            switch version
                case 'case1'
                    state1    = [ 1 log(K_simulA(t-1)/par.K_ss) log(C_simulA(t-1)/par.C_ss)     ];
                    state2    =   dzeta(t)*pi/2 ;
                case 'case2'
                    state1    = [ 1 log(K_simulA(t-1)/par.K_ss) log(C_simulA(t-1)/par.C_ss)     ];
                    state2    =   dzeta(t)*pi/2 ;
                case 'case3'
                    state1    = [ 1 log(K_simulA(t)  /par.K_ss) log(C_simulA(t-1)/par.C_ss)     ];
                    state2    =   dzeta(t)*pi/2 ;
                case 'case4'
                    state1    = [ 1 log(K_simulA(t)  /par.K_ss) log(C_simulA(t-1)/par.C_ss) ...
                                   (log(K_simulA(t)  /par.K_ss))^2                          ...
                                   (log(C_simulA(t-1)/par.C_ss))^2                          ...
                                    log(C_simulA(t-1)/par.C_ss)*log(K_simulA(t)  /par.K_ss)    ];
                    state2    = [ dzeta(t)*pi/2 dzeta(t)^2 dzeta(t)*log(K_simulA(t)/par.K_ss)];
            end
            
            mu_t      = exp(state1*par.eta1)+par.eta_dzeta*sin(state2*par.eta2);
            C_t       = (mu_t)^(-1/par.nu);
            L_t       = (par.b*mu_t*par.labda1*K_simulA(t)^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
            K_next    = par.labda1*K_simulA(t)^par.alpha*L_t^par.beta+(1-par.delta)*K_simulA(t)-C_t;
            % the *_t variables are only used to calculate the conditional expectation

            par.C_t     = C_t;
            par.K_t     = K_simulA(t);
            par.K_next  = K_next;
            
            switch version
                case 'case1'
                    CondExp     = numi_case1(@RHSrealization_sunspot1,par);
                case 'case2'
                    CondExp     = numi(@RHSrealization_sunspot2,par);
                case 'case3'
                    CondExp     = numi(@RHSrealization_sunspot3,par);
                case 'case4'
                    CondExp     = numi(@RHSrealization_sunspot4,par);
            end
            
            C_simulA(t) = CondExp^(-1/par.nu);

            L_simulA(t) = (par.b*C_simulA(t)^(-par.nu)*par.labda1*K_simulA(t)^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
            K_simulA(t+1) = par.labda1*K_simulA(t)^par.alpha*L_simulA(t)^par.beta+(1-par.delta)*K_simulA(t)-C_simulA(t);
            %disp([C_simul(t) C_simulA(t)])
            %pause
        end
        
        figure(3)
        plot([K_simul(1001:1500) K_simulA(1001:1500)])
        figure(4)
        plot([C_simul(1001:1500) C_simulA(1001:1500)])
        mean_Kproj = 100*mean(abs(log(K_simul./K_simulA(1:end))));
        mean_Cproj = 100*mean(abs(log(C_simul./C_simulA(1:end))));
         max_Kproj = 100* max(abs(log(K_simul./K_simulA(1:end))));
         max_Cproj = 100* max(abs(log(C_simul./C_simulA(1:end))));
         disp(' ')
         xxx1 = sprintf('  mean perc. error for K & C:    %7.3f %7.3f',mean_Kproj,mean_Cproj);
         xxx2 = sprintf('   max perc. error for K & C:    %7.3f %7.3f', max_Kproj, max_Cproj);
         xxx0 = sprintf(' standard deviation of K & C:    %7.3f %7.3f',std(log(K_simul)),std(log(C_simul)));
         disp(xxx0)
         disp(xxx1)
         disp(xxx2)
        
end