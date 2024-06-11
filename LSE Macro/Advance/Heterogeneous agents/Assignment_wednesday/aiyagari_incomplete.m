clear;

% Parameters

% Declar parameters to be global such that they can be used in the external
% functions provided.
global a_grid tau mu T alpha delta beta sigma phi N

sigma = 2;              % Coefficient of relative risk-aversion.
beta = 1.05^(-1/4);     % Discount factor.
epsilon = 1e-6;         % Just a small number
N = 200;                % Number of nodes in the grid for wealth.
alpha = 0.3;            % Capital share of income.
delta = 0.025;          % Depreciation rate
mu = 0.5;               % Unemployment benefits

% Grids

phi = 0;            % The debt limit.

% a is individual capital holdings.
a_max = 400;        
a_min = phi;
% The grid is logarithmically spaced in order to put a lot of gridpoints
% where the problem is likely to be highly nonlinear.
a_grid = (exp(linspace(log(a_min+1-a_min),log(a_max+1-a_min),N))-1+a_min)'; % Grid for wealth (endogenous state)

% Transition matrix

T = [0.95, 1-0.95; 0.8, 1-0.8];

% Find the steady state unemployment rate as the eigenvector associated
% with a unit eigenvalue normalized to sum to zero.

[a,b] = eigs(T',1);

a=a/sum(a);

u = a(2);

tau = u/(1-u)*mu;   % Tax rate to balance budget: u*mu*w=(1-u)*tau*w.

% Initial policy function. Pretty stupid guess but it works.

ap_g = a_grid*0.6;

ap_b = a_grid*0.4;

ap_g1 = a_grid*0.6;

ap_b1 = a_grid*0.4;

% Initial bracket for interest rate. This used for the bisection method.

r_bracket = [1/beta-1-epsilon,1/beta-1-epsilon*100];

% A metric to keep check of convergence. Will be used in a while-loop of
% the type while metric_r>1e-6

metric_r = 1;

% Setting up an "almost" Newton-Raphson solver

syms a ap appg appb r w rp wp

euler_g = ((1+r)*a+w*(1-tau)-ap).^(-sigma)-beta*(1+rp)*(T(1,1)*((1+rp)*ap+wp*(1-tau)-appg).^(-sigma)+T(1,2)*((1+rp)*ap+wp*mu-appb).^(-sigma));

euler_b = ((1+r)*a+w*mu-ap).^(-sigma)-beta*(1+rp)*(T(2,1)*((1+rp)*ap+wp*(1-tau)-appg).^(-sigma)+T(2,2)*((1+rp)*ap+wp*mu-appb).^(-sigma));

d_euler_g = diff(euler_g,ap);

d_euler_b = diff(euler_b,ap);

newton_g = ap-euler_g./d_euler_g;

newton_b = ap-euler_b./d_euler_b;

newt_g = matlabFunction(newton_g,'vars',[ a, ap, appg, appb, r, w, rp, wp]);

newt_b = matlabFunction(newton_b,'vars',[ a, ap, appg, appb, r, w, rp, wp]);
    
% r is the midpoint of the bracket.

r = (r_bracket(1)+r_bracket(2))/2;

% w is then implied.

w = (1-alpha)*(((r+delta)/alpha )^(1/(alpha-1)))^(alpha);

% Given w and r we will solve the consumer problem using the solver.

% Set a metric to check convergence

metric = 1;

% Let's rock'n'roll:

while metric>1e-6

    [ap_g1,ap_b1] = newton(r,r,0,0,ap_g1,ap_b1,ap_g,ap_b,newt_g,newt_b);

    metric = max(max(abs([ap_g1-ap_g ap_b1-ap_b])));

    ap_g = ap_g1;

    ap_b = ap_b1;

end

% Done! Consumer problem solved.