function [v1, v2, v3, x, v12, v21, v33] = saveBelmaxp(c,fspace,s)

global beta P r

solvep;

v1  = valfunc1p(c,fspace,s,x);
v2  = valfunc2p(c,fspace,s);
v3  = valfunc3p(c,fspace,s,x);

if nargout > 4
  
v12 =  beta*funbas(fspace, [x, s(:,2)]);

v21 = 0;
v33 = 0;

K = size(P,2);

   for k = 1:K
             
     w  = P(s(:,2),k);
       
     v21  = v21 + bsxfun(@times, w, funbas(fspace, gridmake(s(:,1), k)));       

     v33  = v33 + 1/(1+r)*bsxfun(@times, w, funbas(fspace, gridmake(x, k)));       

     
   end

end
   
