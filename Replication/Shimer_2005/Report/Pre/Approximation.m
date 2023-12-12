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
