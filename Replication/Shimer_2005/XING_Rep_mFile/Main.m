%%%%%% This solves with slight change of Shimer 2005 AER

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
mu_logp = 1;                        % Labor productivity
sigma_logp = 0.05;
rho = 0.8;                          % Persistence
mu_epsilon = 0;                     % Stochastic term of labor productivity
sigma_epsilon = 0.03;
n = 250;                            % Grid size of productivity
[grid_logp, P] = discretizeAR1_Tauchen(0,rho,sigma_epsilon,n,35);
p = exp(grid_logp);                 % p's grid

%% 1) Discretization method 
% Initialize
[U_d, W_d, J_d] = deal(ones(max_iter,n)*tol);           %Value function
[fstar_d, qstar_d, thetastar_d, wstar_d, ustar_d] = deal(ones(n,1)); % Optimal control
[f, q, theta, w] = deal(ones(n,1));                 % Working control
[diffU, diffW, diffJ] = deal(ones(n,1));                % Convergence check

for I = 1:max_iter
    % iterate over productivity grids
    for i = 1:n
        % Working maximized control
        q(i) = c / (delta_ * expec(J_d,P,I,i));
        f(i) = m^(1/alpha) * q(i)^(1 - 1/alpha);
        theta(i) = f(i)/q(i);
        w(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * theta(i);
       
        % Value function iteration
        U_d(I+1,i) = z + delta_ * (f(i) * expec(W_d,P,I,i) + ...
                            (1-f(i)) * expec(U_d,P,I,i));
        W_d(I+1,i) = w(i) + delta_ * ((1-s) * expec(W_d,P,I,i) + ...
                                s * expec(U_d,P,I,i));
        J_d(I+1,i) = p(i) - w(i) + delta_ * (1-s) * expec(J_d,P,I,i);

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
        qstar_d(i) = c / (delta_ * expec(J_d,P,I+1,i));
        thetastar_d(i) = (qstar_d(i) / m)^(-1/alpha);
        fstar_d(i) = m^(1/alpha) * qstar_d(i)^(1 - 1/alpha);
        wstar_d(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * thetastar_d(i);
        ustar_d(i) = delta_ / (delta_ + fstar_d(i));
    end
end

figure
hold on
plot(p,fstar_d);
plot(p,wstar_d);
plot(p,ustar_d);
xlabel("Productivity");
legend("Job finding rate", "Wage", "Unemployment rate");
title('Job Finding Rate, Wage and Unemployment Rate on Productivity by Discretization');
hold off
saveas(gcf,"Discretization.jpg");

%% 2) Parametric Approximation
% Initialize
[U_a, W_a, J_a] = deal(ones(max_iter,n)*tol);
[fstar_a, qstar_a, thetastar_a, wstar_a, ustar_a] = deal(ones(n,1)); % Optimal control
[f, q, theta, w] = deal(ones(n,1));                 % Working control
[diffU, diffW, diffJ] = deal(ones(n,1));            % Convergence check

% We try Hermite Interpolation with pchip
ppU_old = pchip(p, U_a(1,:));
ppW_old = pchip(p, W_a(1,:));
ppJ_old = pchip(p, J_a(1,:));

for I = 1:max_iter
    U_old = ppval(ppU_old,p);
    W_old = ppval(ppW_old,p);
    J_old = ppval(ppJ_old,p);
    % coefs_old = [ppU_old.coefs, ppW_old.coefs, ppJ_old.coefs];
    % Maximization
    for i = 1:n
        % Working maximized control
        q(i) = c / (delta_ * expec(J_a,P,I,i));
        f(i) = m^(1/alpha) * q(i)^(1 - 1/alpha);
        theta(i) = f(i)/q(i);
        w(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * theta(i);
       
        % Value function iteration
        U_a(I+1,i) = z + delta_ * (f(i) * expec(W_a,P,I,i) + ...
                            (1-f(i)) * expec(U_a,P,I,i));
        W_a(I+1,i) = w(i) + delta_ * ((1-s) * expec(W_a,P,I,i) + ...
                                s * expec(U_a,P,I,i));
        J_a(I+1,i) = p(i) - w(i) + delta_ * (1-s) * expec(J_a,P,I,i);
        
        % % Convergence
        % diffU(i) = abs(U_d(I+1,i) - U_d(I,i));
        % diffW(i) = abs(W_d(I+1,i) - W_d(I,i));
        % diffJ(i) = abs(J_d(I+1,i) - J_d(I,i));
    end
    % Fit new parameter
    ppU_new = pchip(p, U_a(I+1,:));
    ppW_new = pchip(p, W_a(I+1,:));
    ppJ_new = pchip(p, J_a(I+1,:));
    
    % Convergence
    U_new = ppval(ppU_new,p);
    W_new = ppval(ppW_new,p);
    J_new = ppval(ppJ_new,p);
    diffU = norm(U_new - U_old);
    diffW = norm(W_new - W_old);
    diffJ = norm(J_new - J_old);
    diff_a = max([diffU,diffW,diffJ]);
    if diff_a <= tol
        break;
    else
        ppU_old = ppU_new;
        ppW_old = ppW_new;
        ppJ_old = ppJ_new;
    end
end

if I == max_iter
    error('NotConverge');
else
    % For the converged value function, derive the optimal control given p
    coefs = [ppU_new.coefs, ppW_new.coefs, ppJ_new.coefs];
    for i = 1:n
        qstar_a(i) = c / (delta_ * expec(J_a,P,I+1,i));
        thetastar_a(i) = (qstar_a(i) / m)^(-1/alpha);
        fstar_a(i) = m^(1/alpha) * qstar_a(i)^(1 - 1/alpha);
        wstar_a(i) = beta_ * p(i) + (1-beta_) * z + beta_ * c * thetastar_a(i);
        ustar_a(i) = delta_ / (delta_ + fstar_a(i));
    end
end

% print("For the approximation method, the coeffecients for U, W and J are given by %6.4f", coefs);
figure
hold on
plot(p,fstar_a);
plot(p,wstar_a);
plot(p,ustar_a);
xlabel("Productivity");
legend("Job finding rate", "Wage", "Unemployment rate");
title('Job Finding Rate, Wage and Unemployment Rate on Productivity by Approximation');
hold off
saveas(gcf,"Approximation.jpg");


%% Q(b)
logp_e = [0.4 0.7 1 1.3 1.6];
p_e = exp(logp_e);                             % Productivity to evaluate
q_e = real(pchip(p,qstar_a,p_e));
f_e = real(m^(1/alpha) * q_e.^(1 - 1/alpha));
w_e = beta_ * p_e + (1-beta_) * z + beta_ * c * f_e ./ q_e;
u_e = delta_ ./ (delta_ + f_e);
Table_b = array2table([logp_e' p_e' w_e' u_e' f_e'], ...
    'VariableNames',{'logp','p','Wage','Unemployment', 'Job Arriving'});
fprintf("The corresponding data are:\n");
disp(Table_b);

%% Q(c)
sd_p = std(p);
sd_u_d = std(ustar_d);
sd_f_d = std(fstar_d);
sd_u_a = std(ustar_a);
sd_f_a = std(fstar_a);
Table_c = array2table([sd_p,sd_u_d,sd_f_d,sd_u_a,sd_f_a], ...
    "VariableNames",{'P','Discre std. of Unemploy', 'Discre std. of Job finding','Approxi std. of Unemploy','Approxi std. of Job finding'});
fprintf("The approximation model gives standard deviation\n");
disp(Table_c)