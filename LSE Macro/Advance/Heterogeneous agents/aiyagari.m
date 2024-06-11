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

% Done. 

% Now let's solve the full problem. Outer loop is for the interest rate.

while metric_r>1e-6
    
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

    % Done, that's it. Consumer problem solved. Find the distribution.

    M = transition(ap_g,ap_b);
    
    opts.disp = 0;
    [V,D] = eigs(M',1,1,opts);
    V = V/sum(V);
    Pe = V(1:N,1);
    Pu = V(N+1:end,1);

    D = [Pe,Pu];
    
    % Type plot(a_grid,D) to take a look at it. But you might want to wait
    % until the program has converged.
    
    % Average capital according to the distribution.

    mean_a = sum(D'*a_grid);
    
    % Implied interest rate

    implied_r = alpha*(mean_a/(1-u))^(alpha-1)-delta;
    
    % s is a "smoothing parameter". I don't want to use the bisection
    % method fully, but update the brackets more slowly. If s=1 it's a
    % standard bisection.
    
    s = 0.01;

    if implied_r>r
        r_bracket(2) = s*r+(1-s)*r_bracket(2);
    else
        r_bracket(1) = s*r+(1-s)*r_bracket(1);
    end
    
    % New r.

    r_new = (r_bracket(1)+r_bracket(2))/2;
    
    % Check convergence.

    metric_r = abs(r-implied_r)*400;
    
    % Print the bracket.
    
    [r_bracket*400 implied_r*400 metric_r]
    
    % And finally update r

    r = r_new;

end

% Done! The Aiyagari model is solved.

% Now let's calculate an impulse response

% Time when we're back to the steady state

time_T = 300;

% Persistence of shock

rho = 0.95;

% Initialize vector of TFP.

z = zeros(time_T,1);

z(1) = -0.01;

for i = 1:time_T
    
    z(i+1) = z(i)*rho^i;
    
end

z(end) = 0;

% Initialize a vector of expected capital throughout the transition.
% A pretty stupid guess, but it will do.

A = ones(time_T+1,1)*mean_a;

% Set a metric again.

metric_A = 1;

while metric_A>1e-6
    
    % Calculate interest rates for the vectors z and A:

    R = alpha*exp(z).*(A./(1-u)).^(alpha-1)-delta;
    
    % Create some storage space for the policy functions.

    apg_imp = zeros(N,time_T);

    apb_imp = zeros(N,time_T);
    
    % We know that we're back in the steady state in period time_T.

    apg_imp(:,end) = ap_g;

    apb_imp(:,end) = ap_b;
    
    % Solve the problem backwards:

    for i = 2:time_T

        [apg_imp(:,time_T-i+1),apb_imp(:,time_T-i+1)] = newton(R(time_T-i+2),R(time_T-i+3),z(time_T-i+2),z(time_T-i+3),apg_imp(:,time_T-i+2),apb_imp(:,time_T-i+2),apg_imp(:,time_T-i+2),apb_imp(:,time_T-i+2),newt_g,newt_b);

    end
    
    % Initial distribution:

    dist = V';
    
    % Update it for each time period and calculate the new sequence of A.

    for i = 1:time_T
        M = transition(apg_imp(:,i),apb_imp(:,i));
        dist = dist*M;
        A_new(i,1) = dist*[a_grid;a_grid];
    end
    
    % Check convergence:

    metric_A = max(abs(A(1:end-1)-A_new))
    
    % Update slowly!

    s = 0.1;

    A = s*A_new+(1-s)*A(1:end-1);

    A(end+1) = mean_a;

end

% Done! Now let's plot the results!

plot(100*(log([mean_a;A])-log(mean_a)))

Capital = [mean_a;A(1:end-1)];
Kss = mean_a;
Output = exp(z).*Capital.^(alpha).*(1-u).^(1-alpha);
Yss = exp(0).*Kss.^(alpha).*(1-u).^(1-alpha);
Investment = Capital(2:end)-(1-delta)*Capital(1:end-1);
Iss = mean_a-(1-delta)*mean_a;
Consumption = Output(1:end-1)-Investment;
Css = Yss-Iss;

subplot(2,2,1);
plot(100*(log([Yss;Output(1:20)])-log(Yss)),'LineWidth',1.5,'color','k');
title('Output, $Y_t$','Interpreter','latex')
axis tight
subplot(2,2,2);
plot(100*(log([Css;Consumption(1:20)])-log(Css)),'LineWidth',1.5,'color','k');
title('Consumption, $C_t$','Interpreter','latex')
axis tight
subplot(2,2,3);
plot(100*(log([Iss;Investment(1:20)])-log(Iss)),'LineWidth',1.5,'color','k');
title('Investment, $I_t$','Interpreter','latex')
axis tight
subplot(2,2,4);
plot(100*([0;z(1:20)]),'LineWidth',1.5,'color','k');
title('Productivity, $Z_t$','Interpreter','latex')
axis tight
suptitle('Impulse Response Functions');
    