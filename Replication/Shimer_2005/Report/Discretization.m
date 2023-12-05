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