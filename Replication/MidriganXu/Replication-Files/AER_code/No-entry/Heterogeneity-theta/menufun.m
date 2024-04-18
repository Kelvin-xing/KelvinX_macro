function [out1, out2] = menufun(flag,s,x)

global r smin smax fspace cr egrid eta alpha delta W 

n = size(s,1);
 
A = s(:,1);

switch flag
    
    case 'b';   
        
E    = exp(egrid(s(:,2)));
rr   = funeval(cr,fspace,s);

L    = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K    = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;
Y    = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D    = Y - W.*L - (r+delta).*K; 
    
out1 = smin(1)*ones(n,1);  
out2 = min(D + (1+r)*A - eps^(1/4), smax(1));
    

    
    case 'f';
        
E    = exp(egrid(s(:,2)));
rr   = funeval(cr,fspace,s);

L    = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K    = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;
Y    = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D    = Y - W.*L - (r+delta).*K; 
C    = D + (1+r)*A - x;     

out1 = log(C);


end
