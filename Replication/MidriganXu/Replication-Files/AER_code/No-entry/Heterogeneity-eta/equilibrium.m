function y = equilibrium(x, Kbar, Lbar)

global eta1 eta2 eta3 alpha Perg egrid

W = x(1);
R = x(2); 

% Compute total losses from misallocation

eta = eta1; 

L1   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K1   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y1   = exp(egrid).^(1-eta).*(L1.^alpha.*K1.^(1-alpha)).^eta;

eta = eta2; 

L2   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K2   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y2   = exp(egrid).^(1-eta).*(L2.^alpha.*K2.^(1-alpha)).^eta;

eta = eta3; 

L3   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K3   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y3   = exp(egrid).^(1-eta).*(L3.^alpha.*K3.^(1-alpha)).^eta;

y = zeros(2,1);

y(1) = Lbar - 1/3*Perg'*(L1 + L2 + L3);
y(2) = Kbar - 1/3*Perg'*(K1 + K2 + K3);
