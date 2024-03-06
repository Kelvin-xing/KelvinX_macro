%%**************************************************************************
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
%  Parameterize modified conditional expectation (= to current k choice)
%--------------------------------------------------------------------------

%Cleaning
clear all; close all; clc

%%
%     newpeaproi1.m MATLAB program to solve the neoclassical growth model
%
%     files needed to run this program:
%     1. RHSrealization_K.m (function to calculate thing inside E_t[]
%     2. rssfcn.m           (function to calculate rss for nonlinear regression
%     3. makepoly.m         (function to generate tensor product polynomialterms
%     4. numi.m             (function to do numerical integration
%
%============================================================================

global T1 T po_k po_t dfactor nu alpha rho delta cs lnts ks ks_S ts_mean ts_std sigma psi

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
%     'previous' for initial values from initC (stochastic simul PEA)
%                %   initial values from initD (quadrature simul PEA)
                   
po_k    =  2;           % Order of Polynomial for k
po_t    =  2;           % Order of Polynomial for theta


maxiter = 2000;         % Maximum Iterations to find
            			% parameters of polynomial
psitol  = 1e-6;         % Convergence criterion
lrate   = 0.1;            % parameter to control updating of coefficient
            			% lrate = 1 means complete updating
            			% lrate < 1 means partial updating
			
% 1.3 Compute the Steady State
% ----------------------------

k_ss      = ( (1-dfactor+delta*dfactor) / (alpha*dfactor) )^(1/(alpha-1));
c_ss      = k_ss^alpha-delta*k_ss;


%% 2 Some Preparatory Work
% =======================

% 2.1 Simulate exogenous process
randn('state',20110629) %State of generator
lnts      = zeros(T,1);
epsi      = randn(T,1)*sigma;

lnts(1)   = 0;
for ti  = 2:T
	lnts(ti)= rho*lnts(ti-1) + epsi(ti);
end;
ts = exp(lnts);
ts_mean = mean(ts(T1:T-1));
ts_std  =  std(ts(T1:T-1));

%The normalizezd value actually being used is
ts_S = (ts-ts_mean)/ts_std;

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
        case 'fresh'
            load dataB % use data in dataB to come up with good starting values
                        
            psi = [ ??? ??? ??? ...
                    ??? ??? ??? ...
                    ??? ??? ???    ]';
            ks_mean =  ???;
            ks_std  =  ???;
        case 'previous'
            switch quadrature
                case 'none'
                    load initC %has psi,ks_mean,ks_std,ks_mean_new,ks_std_new
                case 'herm'
                    load initD %has psi,ks_mean,ks_std,ks_mean_new,ks_std_new
            end
    end

    
%% 4 The Main Loop when no numerical integration is used
% =====================================================

for psiter      = 1:maxiter;

    % 4.1 Generate simulated series
    % -----------------------

    ks(1)   = k_ss;

    for ti  = 2:T
        %construct scaled state variables
        ks_S(ti-1) = (ks(ti-1)-ks_mean)/ks_std;
        polynomial = makepoly([po_k po_t],[ks_S(ti-1) ts_S(ti)]);
    	ks(ti)     = polynomial*psi;
        cs(ti)     = ts(ti)*ks(ti-1)^alpha+(1-delta)*ks(ti-1)-ks(ti);
    end;

%   ks_mean_new = mean(ks(T1-1:T-2));
%   ks_std_new  =  std(ks(T1  :T-1));
%
%    I found that the algorithm works much better if the scaling of the
%    variables is kept fixed across iterations. This is the reason why
%    there is no updating of the mean and standard deviation
%    This does imply that the scaling may not be optimal 
%    (unless you run the algorithm a couple times). But remember that the
%    scaling is only used to make the explanatory terms more orthogonal
%
%

    ks_S = (ks-ks_mean)/ks_std;
    
    % 4.2 Run the regression
    % recall that ks(ti) is end-of-period ti capital
    
    switch quadrature
        case 'none'
        %XXX complete the following expression
        Y = ;    
        Y = ks(T1:T-1).*Y;
        case 'herm'
        %XXX in RHSrealization complete expressions
        Y = numi(@RHSrealization_K,N_herm);            
    end

    X = makepoly([po_k po_t],[ks_S(T1-1:T-2) ts_S(T1:T-1)]);

    psin  = (X'*X)\X'*Y;            

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
        ks_mean_new = mean(ks(T1-1:T-2));
        ks_std_new  =  std(ks(T1  :T-1));
        save initC psi ks_mean ks_std ks_mean_new ks_std_new
        psiC = psi; 
    case 'herm'
        ks_mean_new = mean(ks(T1-1:T-2));
        ks_std_new  =  std(ks(T1  :T-1));
        save initD psi ks_mean ks_std ks_mean_new ks_std_new
        psiD = psi; 
end

% The following end statement ends the psiter do-loop.
end;

%**************************************************************************
%                        End of Main Program                          
%**************************************************************************


