function [v, x, vd] = saveBelmaxw(c, fspace, s)

global beta Pw

solvew; 
v  = valfuncw(c, fspace, s, x);

if nargout > 2
  
  vd = 0;
   
  K  = size(Pw, 2);

  for k = 1 : K
              
     w   = Pw(s(:,2),k);
     vd  = vd  + beta*bsxfun(@times, w, funbas(fspace,gridmake(x, k)));      
     
  end
  
end