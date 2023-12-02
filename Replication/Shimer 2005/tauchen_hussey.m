function [grid_lambda, P] = tauchen_hussey(rho, mu_lambda, sigma_lambda, mu_epsilon, sigma_epsilon, n)
    % Given parameters of a Normal distribution AR(1) process, produce a
    % grid with size n, and corresponding transition matrix P according to
    % Tauchen-Hussey 1991's method, following Kopecky(2007)'s lecture note
    grid_lambda = ones(n,1);
    sigma = sigma_epsilon / sqrt(1-rho^2);  
    m = round(mu_lambda / sigma) - 1;
    grid_lambda(1) = mu_lambda - 
end