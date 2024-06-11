%%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
% sunspotlearning_linear1.m
%
% this program is similar to the log-linear sunspot case of 
% McGough, Meng, and Xue (2011) in the following way:
% the way the sunspot enters makes it possible that the solution for 
% consumption is identical to the one of the loglinear case (which can
% be used as initial conditions)
%
%% overview program
%
%  1. set parameter values algorithm
%  2. set parameter values model & steady state
%  3. get the linearized system
%  4. generate time series with linearized system
%  5. get the nonlinear PEA approximation

clc
clear all

%% 
%
%   1. Parameter values algorithm
%

T              = 11000;
T1             = 1000;
conv_criterion = 1E-8;
dampening      = 0.5;

initialvalues        = 'linear'; % use either linear or previous
                
options = optimset('Display','off','MaxFunEvals',1E5,'MaxIter',...
                  1E5,'TolFun',1E-10,'TolX',1E-10);              
% these are parameters for the minimization routine
% this is a numerical problem within the bigger numerical problem. 
% such innerloop problems should be done more accurately
%
%% 
%   2. Parameter values model
%

par.a       = 0.3;              % capital share production        
par.b       = 0.7;              % labor share production
par.delta   = 0.025;            % depreciation
par.dfactor = 0.99;             % discount factor
par.beta    = 1.1;              % aggregate labor coefficient
par.chi     = 0.08;             % disutility of labor
        
% the following two parameters are adjusted in MMX
        
par.alpha   = 0.25;             % aggregate capital coefficient        
par.nu      = 0.1;              % risk aversion (sigma in MMX)

par.vola    = 0.01;             % standard deviation of the shock
        
% labda1 & labda2 are scalingscoefficients to get L_ss = K_ss = 1
% note that these scalingscoefficients only affect the steady state
% and not of the dynamic properties of the model
        
L_ss        = 1;
K_ss        = 1;
par.labda1  = (1-par.dfactor*(1-par.delta))/(par.dfactor*par.a);
Y_ss        = par.labda1;
C_ss        = Y_ss - par.delta*K_ss;
par.labda2  = par.b*par.labda1*C_ss^(-par.nu);
        
% check steady state (these numbers should be 0)   
disp('check steady state')
disp([-C_ss^(-par.nu)+par.dfactor*C_ss^(-par.nu)*(par.a*par.labda1+1-par.delta)   ...
      -par.labda2+par.b*par.labda1*C_ss^(-par.nu)                                            ...
	   par.labda1+(1-par.delta)-C_ss-1                                          ]);
        

        
 %%
 %   3. Get log-linearized system
 %
 
% linearize budget constraint: k(t+1) = dd_c*c(t) + dd_k*k(t) + dd_l*l(t)
dd_c = -C_ss/K_ss;
dd_k = (1-par.delta)+par.alpha*Y_ss/K_ss;
dd_l =              +par.beta *Y_ss/K_ss;

% linearize FOC labor: l(t) = ff_k*k(t) + ff_c*c(t)
ff_k =  par.alpha/(par.chi+1-par.beta);
ff_c = -par.nu   /(par.chi+1-par.beta);

% substitute out labor from budget constraint:
d_c  = dd_c + dd_l*ff_c;
d_k  = dd_k + dd_l*ff_k;

%linearize Euler equation: c(t) = bb_c*c(t+1) + bb_k*k(t+1) + bb_l*l(t+1)
bb_c =  1;
bb_k = -(par.alpha-1)*par.dfactor*par.a*(Y_ss/K_ss)/par.nu;
bb_l =  -par.beta    *par.dfactor*par.a*(Y_ss/K_ss)/par.nu;

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
%    4. Generate time series using linearized solution
%

 x_tilde = zeros(T,2);
 X       = zeros(T,2);

% generate sunspot variable

rng(20110108); % this replace the former Matlab command to set the seed
% if you have an older version of Matlab replace it with 
% randn('seed', 20110108);
dzeta = par.vola*randn(T,1);
 
 
% generate time paths

for t = 2:T
    x_tilde(t,:) = (JJ*x_tilde(t-1,:)' + [0;dzeta(t)])';
end

X = [x_tilde(:,1) + log(K_ss), ... 
     x_tilde(:,2) + log(C_ss)];
X =  exp(X);

%%
%   5. PEA
%

%
% preparation for PEA
%

switch initialvalues
    case 'linear' 
        % use the linearized solution (for log(C(t)) to get initial conditions for the policy
        % function of period t marginal utility:
        eta = [-par.nu*log(C_ss);-par.nu*JJ(2,1);-par.nu*JJ(2,2);-par.nu];
    case 'previous'
        load eta_linear1
end

C_simul = zeros(T  ,1);
L_simul = zeros(T  ,1);
Y_simul = zeros(T  ,1);
K_simul = zeros(T+1,1);
K_temp  = zeros(T+1,1);

K_simul(1) = K_ss;
K_simul(2) = K_ss;
C_simul(1) = C_ss;

 error_iter = 1000;

 % start the PEA iteration loop
 while error_iter > conv_criterion

     % generate time series
     % XYZ complete this part yourself
     
    for t = 2:T;
        %generate C(t), K(t), & L(t):
        
        % define state used to make forecast:
        % (note that MMX use k(t-1) not k(t))
        state        = ??? ;
        C_simul(t)   = ??? ;

%        the following four lines give solutions for K & L according to linearized solution         
%        k_tilde      = d_k *log(K_simul(t-1)/K_ss)+d_c *log(C_simul(t-1)/C_ss);
%        K_temp(t)    = exp(log(K_ss) + k_tilde);
%        l_tilde    = ff_k*log(K_simul(t)  /K_ss)+ff_c*log(C_simul(t)  /C_ss);
%        L_simul(t) = exp(log(L_ss) + l_tilde);

        L_simul(t)   = ??? ;
        K_simul(t+1) = ??? ;
        
        % finally generate the thing inside the conditional expectation:
        Y_simul(t) = ??? ;
    end

    % do the regression step
   
    Y_simul = Y_simul(T1:T); 
    X_simul = [ones(T-T1+1,1) log(K_simul(T1-2:T-2)/K_ss) log(C_simul(T1-2:T-2)/C_ss)];
    X_dzeta =  dzeta(T1-1:T-1); 

    % give initial values for 3 free values of eta and set fixed parameter
    eta_in     = eta(1:3);
    eta_dzeta  = eta(4);
    % inputs of rssfcn are eta,Y_simul,X_simul,X_dzeta, and eta_dzeta and
    % minimization is with respect to eta
    eta_hat    = fminsearch(@(eta) rssfcn(eta,Y_simul,X_simul,X_dzeta,eta_dzeta),eta_in,options);    
    error_iter = sum(abs(eta_hat-eta(1:3)));
    xxx1 = sprintf('error = %4e',error_iter);
    disp(xxx1)
    xxx2 = sprintf('eta   = %7.3f %7.3f %7.3f %7.3f',eta);  
    disp(xxx2)

    %update beliefs using actual law of motion and previous beliefs
    
    eta(1:3)   = (1-dampening)*eta_hat + dampening*eta(1:3);

 end
 
 save eta_linear1 eta