% compute factor choices and output given (A)

A = exp(s(:,1));
E = exp(s(:,2));

a = max(r,-0.05)*ones(ns,1);
b = 25*ones(ns,1);

rr = a;

L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta));
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta));

sa = sign(lambda*A + (lambda-1)*kappa./E - K);   % this must be < 0

rr = b;

L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta));
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta));

sb = sign(lambda*A + (lambda-1)*kappa./E - K);   % this must be > 0

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
      
L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta));
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta));

sx = sign(lambda*A + (lambda-1)*kappa./E - K);   % this must be > 0
      
       x = x - sx.*dx;
     
  end

Y = (L.^alpha.*K.^(1-alpha)).^eta;
D = Y - W.*L - (r+delta).*K;  

cr = funfitxy(fspace,s,log(rr));
cd = funfitxy(fspace,s,log(D));

clear x;