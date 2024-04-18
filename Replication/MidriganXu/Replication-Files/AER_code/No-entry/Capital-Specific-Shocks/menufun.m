function [out1, out2] = menufun(flag,s,x)

global r smin smax fspace cr egrid eta alpha delta W gamma

n = size(s,1);
 
A = s(:,1);

switch flag
    
    case 'b';   
        
E    = exp(egrid(s(:,2)));
rr   = funeval(cr,fspace,s);

L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*((rr+delta)./((1-alpha).*E.^(1-1/gamma))).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
Y = eta^(eta/(1-eta))*(alpha^gamma*W^(1-gamma) + ((1-alpha)).^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*eta/(eta-1));  

D    = Y - W.*L - (r+delta).*K; 
    
out1 = smin(1)*ones(n,1);  
out2 = min(D + (1+r)*A - eps^(1/4), smax(1));
    

    
    case 'f';
        
E    = exp(egrid(s(:,2)));
rr   = funeval(cr,fspace,s);

L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
K = eta^(1/(1-eta))*((rr+delta)./((1-alpha).*E.^(1-1/gamma))).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta)));
Y = eta^(eta/(1-eta))*(alpha^gamma*W^(1-gamma) + ((1-alpha)).^gamma.*E.^(gamma-1).*(rr+delta).^(1-gamma)).^(1/(1-gamma)*eta/(eta-1));  

D    = Y - W.*L - (r+delta).*K; 
C    = D + (1+r)*A - x;     

out1 = log(C);


end
