function [v1, v2, x, v12, v21] = saveBelmaxt(c,fspace,s)

global beta P fspaceu cu sminu kappau cinj fspacei

ns = size(s,1);

solvet;

v1  = valfunc1t(c,fspace,s,x);
v2  = valfunc2t(c,fspace,s);

if nargout > 3
  
v12 =  beta*funbas(fspace, [x, s(:,2)]);

v21  = 0;
vtu  = 0;
vtt  = 0;
        
K = size(P,2);
inj = funeval(cinj, fspacei, s);

      for k = 1 : K
        
        w  = P(s(:,2),k);        
        
        gu  = gridmake(s(:,1) + inj - kappau, k);           
        vtu = vtu + w.*funeval(cu(:,1),fspaceu,gu);
        
        gt  = gridmake(s(:,1), k);           
        Phi = funbas(fspace, gt);
        vtt = vtt + w.*(Phi*c(:,1));

        v21 = v21 + bsxfun(@times, w, Phi);  
        
      end

      cond = vtu>=vtt & gu(:,1) > sminu(1);    % can't move unless sufficiently large stock of assets amin

      v21  = v21.*repmat(~cond, 1, ns);

end
