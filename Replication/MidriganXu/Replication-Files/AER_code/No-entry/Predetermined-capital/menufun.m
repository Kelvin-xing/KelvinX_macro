function [out1, out2] = menufun(flag,s,x)

global r smin smax fspace cd

n = size(s,1);
 
A = s(:,1);

switch flag
    
    case 'b';   
        

D    = funeval(cd,fspace,s);
    
out1 = smin(1)*ones(n,1);  
out2 = min(D + (1+r)*A - eps^(1/4), smax(1));
    
    case 'f';
        
D  = funeval(cd,fspace,s);
C  = D + (1+r)*A - x;     

out1 = log(C);


end
