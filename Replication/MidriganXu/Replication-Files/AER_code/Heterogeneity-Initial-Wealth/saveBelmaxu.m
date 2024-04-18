function [v1, v2, v3, x, v12, v21, v33] = saveBelmaxu(c,fspace,s)

global beta P fspacep cp sminp kappap r

ns = size(s,1);

solveu;

v1  = valfunc1u(c,fspace,s,x);
v2  = valfunc2u(c,fspace,s);
v3  = valfunc3u(c,fspace,s,x);


if nargout > 4
  
v12 =  beta*funbas(fspace, [x, s(:,2)]);

v21  = 0;
v33  = 0;

vuu  = 0;
vup  = 0;
        
K = size(P,2);

      for k = 1 : K
        
        w  = P(s(:,2),k);        
        
        gp  = gridmake(s(:,1) - kappap, k);           
        vup = vup + w.*funeval(cp(:,1),fspacep,gp);
        
        gu  = gridmake(s(:,1), k);           
        Phi = funbas(fspace, gu);
        vuu = vuu + w.*(Phi*c(:,1));

        v21 = v21 + bsxfun(@times, w, Phi);  
        
        Phi = funbas(fspace, gridmake(x, k));
       
        v33 = v33 + 1/(1+r)*bsxfun(@times, w, Phi);
        
      end

      cond = vup>=vuu & gp(:,1) > sminp(1);    % can't move unless sufficiently large stock of assets amin

      v21  = v21.*repmat(~cond,1,ns);
      v33  = v33.*repmat(~cond,1,ns);
      
end
