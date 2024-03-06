f%%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part I: The Essentials
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%**************************************************************************
%  Wednesday Assignment
%  New simulations parameterized expectations
%  Parameterize usual conditional expectation (= to current marginal utility)
%--------------------------------------------------------------------------

%Cleaning
clear all; close all; clc

%%
%     newpeaproi0.m MATLAB program to solve the neoclassical growth model
%
%     files needed to run this program:
%     1. RHSrealization.m (function to calculate thing inside E_t[]
%     2. rssfcn.m         (function to calculate rss for nonlinear regression
%     3. makepoly.m       (function to generate tensor product polynomialterms
%     4. numi.m           (function to do numerical integration
%
%============================================================================

%%
% File structure:
%
% part 1: 
% set parameters for the model and the algorithm & calculate steady state
%
% part 2: 
% generate a set of realizations (should be fixed across iterations
% so are calculated outside loop
%
% part 3:
% set up the loop over the two choices (non-quadrature & quadrature) and
% set initial values (fresh or read from file)
%
% part 4:
% the main iteration loop to do PEA
%
% part 5:
% plots consumption and capital calculated using the two procedures
%
% part 6: (is in separate file newpeaproi0B.m)
% uses the solution calculated above and calculates the solution for bond prices
% it calculates them B_iter times and then calculates some summary
% statistics. 
% finally it plots those summary statistics. On the x-axis you have the
% number of the sample; in blue you have the outcome for the statistic for
% the different maturities for the non-quad PEA and in black the
% corresponding numbers for the quad PEA

%%
global T1 T po_k po_t dfactor nu alpha rho delta cs lnts ks ks_S sigma psi maturity iter psi_bond
% maturity, iter, & psi_bond are not used in this file but are used in
% newpeaproi0B.m which should be appended to this one

%% 1. Initialize the parameters
% ============================

% 1.1 Initialize model parameters
% -------------------------------

alpha   = 0.33;         % Capital share 
dfactor = .99;          % Time discount factor
nu      = 3.0;          % Risk aversion parameter
delta   = 0.025;        % Depreciation rate

sigma   = .02;          % Standard Deviation for log noise in technology shock
rho     = 0.95;         % Persistence of log technology shock

T       =  2000;        % Total length of simulation
T1      =   501;        % First observation used


% 1.2 Initialize algorithm parameters
% -----------------------------------

N_herm  = 5;            
% number of Hermite nodes for numerical integration

initials = 'fresh';        
%Type 'fresh' for fixed set of initial values 
%     'previous' for initial values from initA (stochastic simul PEA)
%                %   initial values from initB (quadrature simul PEA)
                   
po_k    =  1;           % Order of Polynomial for k
po_t    =  1;           % Order of Polynomial for theta


maxiter = 2000;         % Maximum Iterations to find
            			% parameters of polynomial
psitol  = 1e-6;         % Convergence criterion
lrate   = 0.7;          % parameter to control updating of coefficient
            			% lrate = 1 means complete updating
            			% lrate < 1 means partial updating
			
options = optimset('Display','off','MaxFunEvals',1E5,'MaxIter',...
                  1E5,'TolFun',1E-10,'TolX',1E-10);

% 1.3 Compute the Steady State
% ----------------------------

k_ss      = ( (1-dfactor+delta*dfactor) / (alpha*dfactor) )^(1/(alpha-1));
c_ss      = k_ss^alpha-delta*k_ss;


%% 2 Some Preparatory Work
% =======================

% 2.1 Simulate exogenous process

% Option 1: If you have a recent version of Matlab use

rng(20110629,'twister')
epsi      = randn(T,1)*sigma;

% Option 2: If you have an older version of Matlab set the seed with randn('seed',666);

load shocks
epsi = epsi*sigma;

%it is important you use this seed since the simulated series generated
%with an accurate solution that you will download below is based on this seed.

lnts      = zeros(T,1);

lnts(1)   = 0;
for ti  = 2:T
	lnts(ti)= rho*lnts(ti-1) + epsi(ti);
end;
ts = exp(lnts);
ts_S = lnts;

% 2.2 Initialize vectors
% -----------------------
cs      = zeros(T,1);
ks      = zeros(T,1);
ks_S    = zeros(T,1);


%**************************************************************************
%                        M a i n   P r o g r a m                          
%**************************************************************************

%% 3 Loop over the two possible PEA versions
% =====================================================

for iii = 1:2

% indicate which PEA version is being used
% ----------------------------------------

    if iii == 1
        quadrature = 'none';    
        %'none' for old-fashioned stochastic PEA
    else
        quadrature = 'herm';
        %'herm' for calculating integral
    end

% initialize psi
% --------------------
% we use tensor product polynomials so dimension of psi = (po_k+1)*(po_t+1)

    switch initials
        case 'fresh' % First-order solution of non-quad when nu = 1
                     % for higher-order you cannot simply add zeros since
                     % location of poly terms shift
             psi = [  0.933756226933222   0.284693110705828  -0.530763511673379  -0.178833679612226]';
        case 'previous'
            switch quadrature
                case 'none'
                    load initA   % this contains psi
                case 'herm'
                    load initB   % this contains psi
            end
    end

    
%% 4 The Main Iteration Loop 
% ==========================

for psiter      = 1:maxiter;

    % 4.1 Generate simulated series
    % -----------------------

    ks(1)   = k_ss;
    
    for ti  = 2:T
        %construct scaled state variables
        ks_S(ti-1) = log(ks(ti-1));
        polynomial  = makepoly([po_k po_t],[ks_S(ti-1) ts_S(ti)]);
    	cs(ti)     = exp(-polynomial*psi/nu);
        ks(ti)     = ts(ti)*ks(ti-1)^alpha+(1-delta)*ks(ti-1)-cs(ti);
    end;

    
    % 4.2 Run the regression
    % recall that ks(ti) is end-of-period ti capital
    
    switch quadrature
        case 'none'
        %XXX complete the following expression    
        Y =     
        case 'herm'
        %XXX in RHSrealization complete expressions
        Y = numi(@RHSrealization,N_herm);            
    end

    X = makepoly([po_k po_t],[ks_S(T1-1:T-2) ts_S(T1:T-1)]);

    switch quadrature
        case 'none'
            psin  = fminsearch(@(coef) rssfcn(coef,Y,X),psi,options);    
        case 'herm'
            psin = (X'*X)\X'*log(Y);            
    end
    % note that you are NOT allowed to take the log if you don't
    % numerically integrate RHS and thus still have an additive prediction
    % error

  	delpsi  = norm(psi-psin);
	s1       = sprintf('%3.0d; after %3.0d iterations the Difference in psi was %6.4d.',iii,psiter,delpsi);
    s2       = sprintf('     new estimate for first four elements of psi %f %f %f %f',psi(1:4));
	disp(s1)
    disp(s2)

    if delpsi < psitol
		break;
    else
		psi     = lrate*psin+(1-lrate)*psi;
    end
    
% The following end statement ends the psiter do-loop.
end;

switch quadrature
    case 'none'
        save initA psi 
        psiA = psi;
        csA  = cs;
        ksA  = ks;
        tsA  = ts;
        save dataA csA ksA tsA
    case 'herm'
        save initB psi 
        psiB = psi; 
        csB  = cs;
        ksB  = ks;
        tsB  = ts;
        save dataB csB ksB tsB
end

% The following end statement ends the psiter do-loop.
end;

%**************************************************************************
%                        End of Main Program                          
%**************************************************************************
stop
%% 5 Comparison
% =====================================================

load accurate % contains an accurate solution

figure(1)
plot([ksA(301:600) ksB(301:600) ksAccurate(301:600)])
title('capital')
box off
axis([1 300 18 37])
figure(2)
plot([csA(301:600) csB(301:600) csAccurate(301:600)])
title('consumption')
box off
axis([1 300 1.9 2.65])




