function [out1, out2] = menufunw(flag,s,x)

global r smin smax W lbar 

n = size(s,1);
 
A = s(:,1);

Dw   = W*lbar*(s(:,2) - 1);     

switch flag
    
    case 'bw';   
    
out1 = smin(1)*ones(n,1);  
out2 = min(Dw + (1+r)*A - eps^(1/4),smax(1));

   case 'fw';
    
C    = Dw + (1+r)*A - x;   
out1 = log(C);

end
