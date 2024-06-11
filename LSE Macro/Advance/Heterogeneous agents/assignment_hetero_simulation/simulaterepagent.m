%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%

% this program calculates the decision rule for a rep agent economy and
% simulates a cross-section of agents that all start with a different
% capital stock 
%
% since they all face the same shock the distribution eventually collapses

% with shocks the randomness in z means that even after many periods the
% mass still fluctuates a bit in bins that are close to each other

dynare modelcloglinear.mod   % this generates the policy rule
load dynarerocks             % this uploads the decision rule into the matrix decision
iiicap  = 1;                 % indicates in which column the decision rule for capital is located

logk_ss = decision(1,iiicap);

% construct the grid

N    = 2001;
step = 0.8/(N-1);
gridmin = logk_ss-0.4;
gridmax = logk_ss+0.4;
grid = (gridmin:step:gridmax)';

%
% draw the shocks and generate productivity levels
%
T = 10000;

shock = randn(T,1);

z     = zeros(T,1);
for t = 2:T
    z(t)=rho*z(t-1)+sig*shock(t);
end

% initialize the distribution

hist = ones(N,1)*1/N; % this distributes the total mass (=1) evenly across all grid points

% iterate on the distribution

for t = 1:T
    
    histnew = zeros(N,1);
    
    
%-------------------------------------------------------------------------
%
% Complete below the recursion that specifies the histogram in period.

    for n = 1        :N
        
 
    end
    bar(grid,histnew)
    pause
    hist = histnew;
end
    
