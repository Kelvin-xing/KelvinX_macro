%%%%%% This solves an adapted version of Shimer 2005 AER

%% Initialize
clear;
close;
clc
seed = 1;
rng(seed);
tol = 1e-4;
max_iter = 1100;

%% Q(a)
% Calibration
s = 0.1;                            % Separation rate
r = 0.012;                          % Discount rate
delta_ = 1 / (1 + r);
z = 0.4;                            % Leisure value           
m = 1.355;                          % Average of job flow arrival rate
alpha = 0.72;
beta_ = 0.72;                       % Bargaining power
c = 0.213;                          % Cost of vacancy

% AR(1) process and transition matrix
mu_logp = 0;                        % Labor productivity
sigma_logp = 0.05;
rho = 0.8;                          % Persistence
mu_epsilon = 0;                     % Stochastic term of labor productivity
sigma_epsilon = 0.03;
n = 5;                              % Grid size of productivity
[grid_logp, P] = discretizeAR1_Tauchen(0,rho,sigma_epsilon,n,3);
p = exp(grid_logp);                 % p's grid

%% 1) Discretization method 
% Initialize
[U_d, W_d, J_d] = deal(ones(max_iter,n)*tol);           %Value function
[fstar_d, qstar_d, thetastar_d, wstar_d, ustar_d, vstar_d] = deal(ones(n,1)); % Optimal control
[f, q, theta, w] = deal(ones(n,1));                 % Working control
[diffU, diffW, diffJ] = deal(ones(n,1));                % Convergence check

for I = 1:max_iter
    % iterate over productivity grids
    for i = 1:n
        % Working maximized control
        q(i) = c / (delta_ * J_d(I,:)*P(i,:)');
        f(i) = m^(1/alpha) * q(i)^(1 - 1/alpha);
        theta(i) = f(i)/q(i);
        w(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * theta(i);
       
        % Value function iteration
        U_d(I+1,i) = z + delta_ * (f(i) * W_d(I,:)*P(i,:)' + ...
                            (1-f(i)) * U_d(I,:)*P(i,:)');
        W_d(I+1,i) = w(i) + delta_ * ((1-s) * W_d(I,:)*P(i,:)' + ...
                                s * U_d(I,:)*P(i,:)');
        J_d(I+1,i) = p(i) - w(i) + delta_ * (1-s) * J_d(I,:)*P(i,:)';

        % Convergence
        diffU(i) = abs(U_d(I+1,i) - U_d(I,i));
        diffW(i) = abs(W_d(I+1,i) - W_d(I,i));
        diffJ(i) = abs(J_d(I+1,i) - J_d(I,i));
    end

    % Check convergence
    diff_U = max(diffU);
    diff_W = max(diffW);
    diff_J = max(diffJ);
    diff = max([diff_J, diff_W, diff_U]);
    if diff < tol
        break;
    end
end

if I == max_iter
    error('NotConverge');
else
    % For the converged value function, derive the optimal control given p
    for i = 1:n
        qstar_d(i) = c / (delta_ * J_d(I+1,:)*P(i,:)');
        thetastar_d(i) = (qstar_d(i) / m)^(-1/alpha);
        fstar_d(i) = m^(1/alpha) * qstar_d(i)^(1 - 1/alpha);
        wstar_d(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * thetastar_d(i);
        ustar_d(i) = s / (s + fstar_d(i));
        vstar_d(i) = ustar_d(i) * thetastar_d(i);
    end
end

% Plot Value Function
figure
hold on
plot(p,U_d(I+1,:));
plot(p,W_d(I+1,:));
xlabel("Productivity");
legend("Unemployment value", "Employment value", 'Location','northwest');
title('Value Function by Discretization U&W');
hold off
saveas(gcf,"../Report/Final_report/ValueUWDiscretization.jpg");

figure
plot(p,J_d(I+1,:));
legend("Hiring value", 'Location','northwest');
title('Value Function by Discretization J');
saveas(gcf,"../Report/Final_report/ValueJDiscretization.jpg");

% Plot Optimal Control
figure
hold on
plot(p,fstar_d);
plot(p,wstar_d);
plot(p,ustar_d);
xlabel("Productivity");
legend("Job finding rate", "Wage", "Unemployment rate");
title('Job Finding Rate, Wage and Unemployment Rate on Productivity by Discretization');
hold off
saveas(gcf,"../Report/Final_report/ControlDiscretization.jpg");


%% 2) Parametric Approximation
% Initialize
[U_coef_old, W_coef_old] = deal(ones(1,n)*70);     
J_coef_old = ones(1,n)*tol;                         % Initial guess of value function coefficients
[U_a_old, W_a_old, J_a_old] = deal(ones(1,n));      % Working value function
[U_a_new, W_a_new, J_a_new] = deal(ones(1,n));      % Working value function
[fstar_a, qstar_a, thetastar_a, wstar_a, ustar_a] = deal(ones(n,1)); % Optimal control
[f, q, theta, w] = deal(ones(n,1));                 % Working control
% We do a simple Lagrange Interpolation of the value functions
F_hat = @(p,coef) coef(1) + p * (coef(2) + p * (coef(3) + p * (coef(4) + p * coef(5))));

for I = 1:max_iter
    % Initial guess
    for i = 1:n
        % Working value function
        U_a_old(i) = F_hat(p(i),U_coef_old);
        W_a_old(i) = F_hat(p(i),W_coef_old);
        J_a_old(i) = F_hat(p(i),J_coef_old);
    end
    % Maximization
    for i = 1:n
        % Working maximized control
        q(i) = c / (delta_ * J_a_old*P(i,:)');
        f(i) = m^(1/alpha) * q(i)^(1 - 1/alpha);
        theta(i) = f(i)/q(i);
        w(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * theta(i);
       
        % Value function iteration
        U_a_new(i) = z + delta_ * (f(i) * W_a_old*P(i,:)' + ...
                            (1-f(i)) * U_a_old*P(i,:)');
        W_a_new(i) = w(i) + delta_ * ((1-s) * W_a_old*P(i,:)' + ...
                                s * U_a_old*P(i,:)');
        J_a_new(i) = p(i) - w(i) + delta_ * (1-s) * J_a_old*P(i,:)';
    end

    % Fitting
    U_coef_new = fminsearch(@(coef) err(coef, U_a_new, p), U_coef_old);
    W_coef_new = fminsearch(@(coef) err(coef, W_a_new, p), W_coef_old);
    J_coef_new = fminsearch(@(coef) err(coef, J_a_new, p), J_coef_old);

    % Convergence check
    for i = 1:n
        % Approximation with new coefficients
        U_a_new(i) = F_hat(p(i),U_coef_new);
        W_a_new(i) = F_hat(p(i),W_coef_new);
        J_a_new(i) = F_hat(p(i),J_coef_new);
    end
    % errors
    diffU = norm(U_a_new - U_a_old);
    diffW = norm(W_a_new - W_a_old);
    diffJ = norm(J_a_new - J_a_old);
    diff_a = max([diffU,diffW,diffJ]);
    if diff_a <= tol
        break;
    else
        U_coef_old = U_coef_new;
        W_coef_old = W_coef_new;
        J_coef_old = J_coef_new;
    end
end

if I == max_iter
    error('NotConverge');
else
    % For the converged value function, derive the optimal control given p
    for i = 1:n
        qstar_a(i) = c / (delta_ * J_a_new*P(i,:)');
        thetastar_a(i) = (qstar_a(i) / m)^(-1/alpha);
        fstar_a(i) = m^(1/alpha) * qstar_a(i)^(1 - 1/alpha);
        wstar_a(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * thetastar_a(i);
        ustar_a(i) = s / (s + fstar_a(i));
    end
end
% Show coefficients
coefs = array2table([U_coef_new' W_coef_new' J_coef_new'], ...
    "VariableNames",{'U', 'W','J'}, ...
    "RowNames",{'e'; 'd'; 'c'; 'b'; 'a'});
fprintf("For the approximation method, the coeffecients for U, W and J are given by:\n");
disp(coefs);

% Plot Value Function
figure
hold on
plot(p,U_a_new);
plot(p,W_a_new);
xlabel("Productivity");
legend("Unemployment value", "Employment value", 'Location','northwest');
title('Value Function by Discretization U&W');
hold off
saveas(gcf,"../Report/Final_report/ValueUWApproximation.jpg");

figure
plot(p,J_a_new);
legend("Hiring value", 'Location','northwest');
title('Value Function by Discretization J');
saveas(gcf,"../Report/Final_report/ValueJApproximation.jpg");

% Plot Optimal Control
figure
hold on
plot(p,fstar_a);
plot(p,wstar_a);
plot(p,ustar_a);
xlabel("Productivity");
legend("Job finding rate", "Wage", "Unemployment rate");
title('Job Finding Rate, Wage and Unemployment Rate on Productivity by Approximation');
hold off
saveas(gcf,"../Report/Final_report/ControlApproximation.jpg");


%% Q(b)
logp_e = [-0.10 -0.05 0 0.05 0.10];
p_e = exp(logp_e);                             % Productivity to evaluate
q_e = spline(p,qstar_d,p_e);
f_e = m^(1/alpha) * q_e.^(1 - 1/alpha);
w_e = beta_ * p_e + (1-beta_) * z + beta_ * c * f_e ./ q_e;
u_e = s ./ (s + f_e);
Table_b = array2table([logp_e' w_e' u_e' f_e'], ...
    'VariableNames',{'logp','Wage','Unemployment', 'Job Arriving'});
fprintf("The corresponding data are:\n");
disp(Table_b);

%% Q(c)
% Standard deviation
sd_p = std(p);
sd_u_d = std(ustar_d);
sd_f_d = std(fstar_d);
sd_u_a = std(ustar_a);
sd_f_a = std(fstar_a);

% Correlation Matrix
variables = [ustar_d,vstar_d,fstar_d,p];
corrmain = corrcoef(variables);
corrmatrix = array2table(corrmain,'RowNames',{'u','v','f','p'}, ...
                          'VariableNames',{'u','v','f','p'});
% Table display
Table_cmain = array2table([sd_p,sd_u_d,sd_f_d,sd_u_a,sd_f_a], ...
    "VariableNames",{'P','Discre std. of Unemploy', 'Discre std. of Job finding','Approxi std. of Unemploy','Approxi std. of Job finding'});
fprintf("The model gives standard deviation\n");
disp(Table_cmain)
fprintf("The correlation matrix reads:\n");
disp(corrmatrix)