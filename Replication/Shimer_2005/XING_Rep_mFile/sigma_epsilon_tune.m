%%%%%% This tunes sigma epsilon relative to result in Main.m
%%%%%% Apply only discretization method for simplification

%% Initialize
clear all;
close;
clc
seed = 1;
rng(seed);
tol = 1e-4;
max_iter = 1100;

%% Calibration
s = 0.1;                            % Separation rate
r = 0.012;                          % Discount rate
delta_ = 1 / (1 + r);
z = 0.4;                            % Leisure value           
m = 1.355;                          % Average of job flow arrival rate
alpha = 0.72;
beta_ = 0.72;                       % Bargaining power
c = 0.213;                          % Cost of vacancy

%% AR(1) process and transition matrix
mu_logp = 0;                        % Labor productivity
sigma_logp = 0.05;
rho = 0.8;                          % Persistence
mu_epsilon = 0;                     % Stochastic term of labor productivity
sigma_epsilon = 0.01:0.01:0.05;     % sigma epsilon sample
n = 5;                              % Grid size of productivity
[~, epslistsize] = size(sigma_epsilon);
[sd_p, sd_u_d, sd_f_d] = deal(ones(1,epslistsize));% Allocate size
for j = 1:epslistsize
    [grid_logp, P] = discretizeAR1_Tauchen(mu_logp,rho,sigma_epsilon(j),n,3);
    p = exp(grid_logp);                 % p's grid
    
    %% Value function iteration
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
    
    % Standard deviation
    sd_p(j) = std(p);
    sd_u_d(j) = std(ustar_d);
    sd_f_d(j) = std(fstar_d);
end
% Correlation Matrix with 0.05 sigma epsilon
variables = [ustar_d,vstar_d,fstar_d,p];
corrsigmaep = corrcoef(variables);
corrmatrix = array2table(corrsigmaep,'RowNames',{'u','v','f','p'}, ...
                          'VariableNames',{'u','v','f','p'});

% Display table
Table_csigmaep = array2table([sigma_epsilon', sd_p',sd_u_d',sd_f_d'], ...
    "VariableNames",{'sigma_epsilon','P','Unemploy', 'Job finding'});
fprintf("The model gives standard deviation\n");
disp(Table_csigmaep)
fprintf("The correlation matrix reads:\n");
disp(corrmatrix)