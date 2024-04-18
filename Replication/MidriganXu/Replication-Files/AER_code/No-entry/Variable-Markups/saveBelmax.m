function [v1, v2, x, v12, v21] = saveBelmax(c,fspace,s)

global beta P

ns = size(s,1);

solve;

v1  = valfunc1(c,fspace,s,x);
v2  = valfunc2(c,fspace,s);

if nargout > 3
  
v12 =  beta*funbas(fspace, [x, s(:,2), s(:,2)]);

v21 = 0;

K = size(P,2);

   for k = 1:K
             
     w  = P(s(:,2),k);
       
     v21  = v21 + bsxfun(@times, w, funbas(fspace, [s(:,1), k*ones(ns,1), s(:,2)]));       
     
   end

end
   
