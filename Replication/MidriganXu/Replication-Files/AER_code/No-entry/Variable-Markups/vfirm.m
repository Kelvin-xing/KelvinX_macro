function pai = vfirm(ybar, s)

global egrid P alpha eta r delta W theta

% Compute, given a vector of ybar (predetermined output), expected profits

Exp = zeros(size(s,1), 1); 

for ii = 1 : size(egrid,1)
  
    w  = P(s(:,3),ii); 
  
    Kguess = (alpha/W)^(-alpha)*((1-alpha)/(r+delta))^alpha.*exp(egrid(ii)).^(1-1/eta).*ybar.^(1/eta); % frictionless
    
    mu = W*(1-theta)^(1/alpha).*s(:,1).^(-1/alpha).*(1-alpha)./alpha.*exp(egrid(ii)).^(1/alpha*(1-1/eta)).*ybar.^(1/alpha/eta) - r - delta;
    
    mu(Kguess <= 1/(1-theta)*s(:,1)) = 0;   % unconstrained
    
    Exp = Exp + w.*((alpha + (1-alpha)*(r + delta)./(r + delta + mu)).*(r + delta + mu).^(1-alpha).*exp(egrid(ii)).^(1-1/eta));
    
end


pai = ybar - alpha^(-alpha)*(1-alpha)^(-(1-alpha))*W^alpha.*ybar.^(1/eta).*Exp;

