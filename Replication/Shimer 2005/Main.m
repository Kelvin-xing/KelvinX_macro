%% Initialize
clear;
close;
clc
=seed = 1;
rng(seed);
tol = 1e-6;
max_iter = 1000;

%% Q(a)
% Calibration
s = 0.1;                            % Separation rate
r = 0.012;                          % Discount rate
delta_ = 1 / (1 + r);
z = 0.4;                            % Leisure value           
m = 1.355;                          % Average of job flow arrival rate
alpha = 0.72;
q = @(theta_) m * theta_^(-alpha);  % recruiting rate for firm
f = @(theta_) m * theta_^(1-alpha); % job arrival rate for the unemployed
beta_ = 0.72;                       % Bargaining power
c = 0.213;                          % Cost of vacancy

% AR(1) process and transition matrix
mu_logp = 1;                        % Labor productivity
sigma_logp = 0.015;
rho = 0.8;
mu_epsilon = 0;                     % Stochastic term of labor productivity
sigma_epsilon = 0.03;
n = 100;                           % Grid size of productivity
[grid_logp, P] = tauchen_hussey(rho, mu_logp, sigma_logp, mu_epsilon, sigma_epsilon, n);
p = exp(grid_logp);                 % Each p's value

%% 1) Discretization method 
% Initialize
[U, W, J] = deal(ones(n));
[qstar, thetastar, wstar] = ones(n,1); % Optimal control
[f, q, theta, w] = ones(n,1);          % Working control

% Algorithm: 
% 1) define a grid of p make initial guesses of UWJV of the same size
% 2) for each p, calculate conduct value function iteration of each value function
% 3) fing the optimal q(theta) and w for each p
% 4) calculate the unemployment rate

for I = 1:max_iter 
    % iterate over productivity grids
    for i = 1:n
        % Working control
        q(i) = c / (delta_ * expec(J,P,I,i));
        f(i) = m^(1/alpha) * q(i)^(1 - 1/alpha);
        theta(i) = (q(i) / m)^(-1 / alpha);
        w(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * theta(i);
       
        % Value function iteration
        U(I+1,i) = z + delta_ * (f(i) * expec(W,P,I,i) + ...
                            (1-f(i) * expec(U,P,I,i)));
        W(I+1,i) = w(i) + delta_ * ((1-s) * expec(W,P,I,i) + ...
                                s * expec(U,P,I,i));
        J(I+1,i) = p(i) - w(i) + delta_ * (1-s) * expec(J,P,I,i);
    
        % Check convergence
        diffU = norm(U(I+1,:) - U(I,:));
        diffW = norm(W(I+1,:) - W(I,:));
        diffJ = norm(J(I+1,:) - J(I,:));
        diff = max([diffJ, diffW, diffU]);
        if diff < tol
            break
        end
    end
    if I == max_iter
        raise NotConverge 
    end
    % For the converged value function, derive the optimal control
    qstar(i) = c / (delta_ * expec(J,P,I+1,i));
    thetastar(i) = exp(qstar(i) / m, -1 / alpha);
    wstar(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * thetastar(i);
end

%%% Approximation

%% Q(b)
% logp = [0.4 0.7 1 1.3 1.6];
% p = exp(logp);

%% Q(c)

