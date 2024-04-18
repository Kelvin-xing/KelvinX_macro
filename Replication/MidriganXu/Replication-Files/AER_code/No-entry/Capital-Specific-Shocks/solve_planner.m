function y = solve_planner(x, Ka, La, e)

global eta gamma alpha

L = eta^(1/(1-eta))*(x(1)/alpha)^(-gamma)*(alpha^gamma*x(1)^(1-gamma) + (1-alpha)^gamma.*exp(e).^(gamma-1).*x(2).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*(x(2)./((1-alpha).*exp(e).^(1-1/gamma))).^(-gamma).*(alpha^gamma*x(1)^(1-gamma) + (1-alpha)^gamma.*exp(e).^(gamma-1).*x(2).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));


y = zeros(2,1); 

y(1) = mean(L(:)) - La;
y(2) = mean(K(:)) - Ka;
