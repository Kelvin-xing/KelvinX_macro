  function [v, vup, vuu] = valfunc2u(c,fspace,s)

  global P fspacep cp sminp kappap
 
     vup  = 0;
     vuu  = 0;
        
K = size(P,2);

      for k = 1 : K
        
        w  = P(s(:,2),k);

        gu  = gridmake(s(:,1), k);           
        vuu = vuu + w.*funeval(c(:,1), fspace, gu);

        gp  = gridmake(s(:,1) - kappap, k);           
        vup = vup + w.*funeval(cp(:,1), fspacep, gp);
                      
      end

     cond = vup>=vuu & gp(:,1) > sminp(1);    % can't move unless sufficiently large stock of assets -- Amin
     
     v = vup.*cond + vuu.*(~cond);
     
vup(s(:,1) < kappap + sminp(1)) = vuu(s(:,1) < kappap + sminp(1)) - 5;