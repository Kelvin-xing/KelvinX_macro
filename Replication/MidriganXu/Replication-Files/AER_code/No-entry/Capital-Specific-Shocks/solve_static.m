% compute factor choices and output given (A)

A = s(:,1);
E = exp(egrid(s(:,2)));

a = max(r,-0.05)*ones(ns,1);
b = 25*ones(ns,1);

rr = a;

L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*((rr+delta)./((1-alpha).*E.^(1-1/gamma))).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));

sa = sign(1/(1-theta)*A - K);   % this must be < 0

rr = b;

L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*((rr+delta)./((1-alpha).*E.^(1-1/gamma))).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));

sb = sign(1/(1-theta)*A - K);   % this must be > 0

if any(sa>0)   % unconstrained, don't search further
   
    b(sa>0) = a(sa>0);
    
end

if any(sb<0)
    
    a(sb<0) = b(sb<0);
    
end

% Initializations  

tol = 1e-11;

 dx = 0.5*(b - a);
 x  = a + dx;
 dx = sb.*dx;

% Iteration loop

  while any(abs(dx)>tol)

      dx = 0.5*dx;
      rr = x;
      
L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*((rr+delta)./((1-alpha).*E.^(1-1/gamma))).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));

sx = sign(1/(1-theta)*A - K);   % this must be > 0
      
       x = x - sx.*dx;
     
  end

Y = eta^(eta/(1-eta))*(alpha^gamma*W^(1-gamma) + ((1-alpha)).^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*eta/(eta-1));  
D = Y - W*L - (r+delta)*K;

cr = funfitxy(fspace, s, rr);
cl = funfitxy(fspace, s, L);
ck = funfitxy(fspace, s, K);
cy = funfitxy(fspace, s, Y);
cd = funfitxy(fspace, s, D);


clear x;