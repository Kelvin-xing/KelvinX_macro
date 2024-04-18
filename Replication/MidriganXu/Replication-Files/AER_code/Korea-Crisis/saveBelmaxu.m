function [v1, v2, x, v12, v21] = saveBelmaxu(c,fspace,s)

global beta P

ns = size(s,1);

solveu;

v1  = valfunc1u(c,fspace,s,x);
v2  = valfunc2u(c,fspace,s);


if nargout > 3
  
v12 =  beta*funbas(fspace, [x, s(:,2)]);

v21  = 0;
        
K = size(P,2);

      for k = 1 : K
        
        w  = P(s(:,2),k);        
        
        g  = gridmake(s(:,1), k);           
        Phi = funbas(fspace, g);
        
        v21 = v21 + bsxfun(@times, w, Phi);  
        
      end

 
      
end
