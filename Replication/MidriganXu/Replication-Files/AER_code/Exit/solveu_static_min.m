% compute factor choices and output given (A)

A = nodeunif(150000, -theta*nuk*(kappau) + eps^(1/4), 10);
E = exp(egrid(1) + phiu);                          % lower bound on E

a = max(r,-0.05)*ones(length(A),1);
b = 25*ones(length(A),1);

rr = a;

L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;

sa = sign(1/(1-theta)*A + theta*nuk/(1-theta)*(kappau) - K);   % this must be < 0

rr = b;

L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;

sb = sign(1/(1-theta)*A + theta*nuk/(1-theta)*(kappau) - K);   % this must be > 0

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

rr = a;   % in case unconstrained so never loop
L  = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K  = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;

  while any(abs(dx)>tol)

      dx = 0.5*dx;
      rr = x;
      
L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;

sx = sign(1/(1-theta)*A + theta*nuk/(1-theta)*(kappau) - K);   % this must be > 0
      
       x = x - sx.*dx;
     
  end

Y = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;

Dm = Y - W.*L - (r+delta).*K - W*F;  
Dt = E.*pait;
D = max(Dm,Dt);

Aumin = min(A((1-theta*xai)*D + r*A > 0));