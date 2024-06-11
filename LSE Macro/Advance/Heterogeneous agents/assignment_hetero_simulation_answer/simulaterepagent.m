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
    t
    histnew = zeros(N,1);
    for n = 1        :N
        logk      = grid(n);
        logk_next = logk_ss + decision(2,iiicap)*(logk-logk_ss) + decision(4,iiicap)*z(t);

        % locate the corresponding grid points
        % the lines below are reprogrammed a bit (a mistake regarding the
        % upperbound is corrected and the structure is now more clear)
        
        n_next    = 1+floor((logk_next-gridmin)/step);

        if logk_next >= gridmin && logk_next < gridmax
        weight    = (grid(n_next+1)-logk_next)/(grid(n_next+1)-grid(n_next));
        end
        if logk_next <  gridmin
            n_next = 1;
            weight = 1;
        end
        if logk_next >= gridmax;
           n_next  = N-1;
           weight = 0;
        end

        histnew(n_next)   = histnew(n_next)  +hist(n)*   weight;
        histnew(n_next+1) = histnew(n_next+1)+hist(n)*(1-weight);
 
    end
    bar(grid,histnew)
    pause
    hist = histnew;
end
    
