% Parameters
s = 0.1;
r = 0.012;
z = 0.4;
mu = 1;
sigma = 0.015;
eta = 0.72;
beta = 0.72;
c = 0.213;

% Grid for productivity levels
p_min = exp(mu - 3*sigma);
p_max = exp(mu + 3*sigma);
n_grid = 100;
dp = (p_max - p_min) / (n_grid - 1);
log_p = log(p_min:dp:p_max);

% Value function iteration
max_iter = 1000;
tol = 1e-6;

% Initial guesses for value functions
U = z * ones(n_grid, 1);
W = zeros(n_grid, 1);
J = zeros(n_grid, 1);
V = zeros(n_grid, 1);

for iter = 1:max_iter
    % Value function iteration for each productivity level
    U_new = zeros(n_grid, 1);
    W_new = zeros(n_grid, 1);
    J_new = zeros(n_grid, 1);
    V_new = zeros(n_grid, 1);
    for i = 1:n_grid
        % Arrival rate of jobs and recruiting rate
        theta = exp(log_p(i) - W);
        f = eta * exp(mu + (eta - 1) * log(theta)) .* (theta <= 1);
        q = eta * exp(mu - eta * log(theta)) .* (theta >= 1) .* (theta <= exp(mu));
        
        % Expected value functions
        EU = f .* W + (1 - f) .* U;
        EW = (1 - s) * W + s * U;
        EJ = J;
        EV = q .* J;
        
        % Update value functions
        U_new(i) = z + r * (f' * EU + (1 - f') * U_new);
        W_new(i) = (U_new(i) + beta * J_new(i)) / (1 + beta);
        J_new(i) = exp(log_p(i)) - W_new(i) + r * (1 - s) * EJ' * ones(n_grid, 1);
        V_new(i) = -c + r * q' * EJ;
    end
    
    % Check for convergence
    if max(abs(U_new - U)) < tol && max(abs(W_new - W)) < tol && max(abs(J_new - J)) < tol && max(abs(V_new - V)) < tol
        break;
    end
    
    % Update value functions
    U = U_new;
    W = W_new;
    J = J_new;
    V = V_new;
end

% Equilibrium wage
w = mean(W);

% Plot value functions
figure;
plot(log_p, U, 'r', log_p, W, 'b', log_p, J, 'g', log_p, V, 'm');
legend('Unemployed', 'Employed', 'Firm hiring', 'Firm posting');
xlabel('Log productivity');
ylabel('Value function');