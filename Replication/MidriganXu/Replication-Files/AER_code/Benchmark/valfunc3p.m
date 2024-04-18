  function v = valfunc3p(c, fspace, s, x)

  global P cr egrid phip r  eta alpha delta W
  
v = 0;

K = size(P,2);
     
  for k = 1 : K
              
     w  = P(s(:,2),k);
     
E    = exp(egrid(k) + phip);               % compute tomorrow's dividends
rr   = funeval(cr,fspace,gridmake(x, k));  % tomorrow's multiplier

L    = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K    = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;
Y    = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D    = Y - W.*L - (r+delta).*K;

     
     v  = v  + 1/(1+r)*w.*(funeval(c(:,3), fspace, gridmake(x, k)) + D);       
     
  end

  
