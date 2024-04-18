% compute factor choices and output given (A)

A = s(:,1);
E = exp(egrid(s(:,2)));

E1 = 0*E; % E_{-1} e^... term in capital decision rule

for ii = 1 : k
  
    w  = P(s(:,3),ii); 
  
    E1 = E1 + w.*exp(egrid(ii)).^(1-1/eta);
    
end

yguess = eta^(eta/(1-eta))*(alpha/W)^(alpha*eta/(1-eta))*((1-alpha)/(r+delta))^((1-alpha)*eta/(1-eta))*E1.^(eta/(eta-1));


a = 0.001*min(yguess);
b = 5*max(yguess);

tol = 1e-11;

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;

d  = b - a; 

x1 = a + alpha1*d; 
x2 = a + alpha2*d; 

f1 = vfirm(x1, s); 
f2 = vfirm(x2, s); 

d = alpha1*alpha2*d;

x1new = x1;
x2new = x2;
f1new = f1;
f2new = f2;

while any((d)>tol)

 f1 = f1new;
 f2 = f2new;
 x1 = x1new;
 x2 = x2new;
    
 d = d*alpha2;
 x2new = x1.*(f2<f1)+(x2+d).*(f2>=f1);
 f2new = f1.*(f2<f1)+vfirm(x2+d, s).*(f2>=f1);
  
 x1new = (x1-d).*(f2<f1)+x2.*(f2>=f1);
 f1new = vfirm(x1-d, s).*(f2<f1)+f2.*(f2>=f1);
  
end

ybar = x2new.*(f2new>=f1new)+x1new.*(f2new<f1new);    % Optimal output given e(-1) and a
 
Kguess = (alpha/W)^(-alpha)*((1-alpha)/(r+delta))^alpha.*exp(egrid(s(:,2))).^(1-1/eta).*ybar.^(1/eta); % frictionless
    
mu = W*(1-theta)^(1/alpha).*s(:,1).^(-1/alpha).*(1-alpha)./alpha.*exp(egrid(s(:,2))).^(1/alpha*(1-1/eta)).*ybar.^(1/alpha/eta) - r - delta;
    
mu(Kguess <= 1/(1-theta)*s(:,1)) = 0;   % unconstrained

L  = (alpha/W)^(1-alpha).*((1-alpha)./(r+delta+mu)).^(-(1-alpha)).*exp(egrid(s(:,2))).^(1-1/eta).*ybar.^(1/eta);
K  = (alpha/W)^(-alpha).*((1-alpha)./(r+delta+mu)).^alpha.*exp(egrid(s(:,2))).^(1-1/eta).*ybar.^(1/eta);
Y  = ybar;
D  = Y - w.*L - (r + delta).*K;
rr = r + mu;

Yf = yguess;

cr = funfitxy(fspace, s, rr);

cl = funfitxy(fspace, s, L);
cy = funfitxy(fspace, s, Y);
cd = funfitxy(fspace, s, D); 
ck = funfitxy(fspace, s, K);

cyf = funfitxy(fspace,s,Yf);