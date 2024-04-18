% compute factor choices and output given (A)

A = s(:,1);
E = exp(egrid(s(:,2)));

E1 = 0*E; % E_{-1} e^... term in capital decision rule

for ii = 1 : k
  
    w  = P(s(:,3),ii); 
  
    E1 = E1 + w.*exp(egrid(ii)).^((1-eta)/(1-alpha*eta));
    
end


a = max(r,-0.05)*ones(ns,1);
b = 25*ones(ns,1);

rr = a;

K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E1.^((1-alpha*eta)/(1-eta));

sa = sign(1/(1-theta)*A - K);   % this must be < 0

rr = b;

K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E1.^((1-alpha*eta)/(1-eta));

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
      
K = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E1.^((1-alpha*eta)/(1-eta));

sx = sign(1/(1-theta)*A - K);   % this must be > 0
      
       x = x - sx.*dx;
     
  end

cr = funfitxy(fspace, s, rr);

L = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E.^((1-eta)./(1-alpha*eta)).*E1.^((1-alpha)*eta/(1-eta));
Y = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D = Y - W.*L - (r+delta).*K; 

cl = funfitxy(fspace, s, L);
cy = funfitxy(fspace, s, Y);
cd = funfitxy(fspace, s, D); 
ck = funfitxy(fspace, s, K);
ce = funfitxy(fspace, s, E1); % expectation of current e given e(-1);

clear x;