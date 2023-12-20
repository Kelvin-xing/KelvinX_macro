%%%%%% This sets separation rate as iid draw from normal distribution

%% Initialize
clear;
close;
clc
seed = 1;
rng(seed);
tol = 1e-4;
max_iter = 2000;

%% Q(a)
% Calibration
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
n = 100;                              % Grid size of productivity
[grid_logp, P] = discretizeAR1_Tauchen(0,rho,sigma_epsilon,n,3);
p = exp(grid_logp);                 % p's grid

% Separation Rate
mu_s = 0.1:0.1:0.8;
sigma_s = 0.1;                      % Alternative include 1-z=0.6; 0.075                                    
N = 100;
[~,O] = size(mu_s);
[sd_u_d,sd_f_d,sd_p, sd_s] = deal(ones(O,1));
for o = 1:O
    s = max(normrnd(mu_s(o),sigma_s,N,1),0);
    
    %% Value function iteration.
    % Initialize
    [U_d, W_d] = deal(ones(max_iter,n)*70);           %Value function
    J_d = ones(max_iter,n)*tol;
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
            tempW = 0;
            tempJ = 0;
            % Monte-Carlo Integration
            for iter = 1:N
                tempW = tempW + w(i) + delta_ * ((1-s(iter)) * W_d(I,:)*P(i,:)' ...
                                                + s(iter) * U_d(I,:)*P(i,:)');
                tempJ = tempJ + p(i) - w(i) + delta_ * (1-s(iter)) * J_d(I,:)*P(i,:)';
            end
            W_d(I+1,i) = tempW / N;
            J_d(I+1,i) = tempJ / N;
    
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
            temp = 0;
            for iter = 1:N
                temp = temp + s(iter) / (s(iter) + fstar_d(i));
            end
            ustar_d(i) = temp/N;
            vstar_d(i) = ustar_d(i)*thetastar_d(i);
        end
    end
    % Standard deviation
    sd_p(o) = std(p);
    sd_s(o) = std(s);
    sd_u_d(o) = std(ustar_d);
    sd_f_d(o) = std(fstar_d);
end
%% Correlation Matrix with 0.8 mu_s
variables = [ustar_d,vstar_d,fstar_d,p,s];
corrstocsep = corrcoef(variables);

corrmatrix = array2table(corrstocsep,'RowNames',{'u','v','f','p','s'}, ...
                          'VariableNames',{'u','v','f','p','s'});

% Table display
Table_cstocsep = array2table([mu_s',sd_u_d,sd_f_d,sd_p, sd_s], ...
    "VariableNames",{'mean_s','Unemploy', 'Job-finding','p','s'});
fprintf("The model gives standard deviation\n");
disp(Table_cstocsep)
fprintf("The correlation matrix reads:\n");
disp(corrmatrix)