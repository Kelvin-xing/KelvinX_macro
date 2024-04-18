  function v = valfunc3u(c, fspace, s, x)

  global P cr egrid phip phiu r  eta alpha delta W cp fspacep kappap sminp
  
v = 0;

K = size(P,2);
       
vup  = 0;
vuu  = 0;
        
K = size(P,2);

      for k = 1 : K
        
        w  = P(s(:,2),k);

        gu  = gridmake(x, k);           
        vuu = vuu + w.*funeval(c(:,1), fspace, gu);

        gp  = gridmake(x - kappap, k);           
        vup = vup + w.*funeval(cp(:,1), fspacep, gp);
                      
      end

     cond = vup>=vuu & gp(:,1) > sminp(1);    % can't move unless sufficiently large stock of assets -- Amin
     
v = 0;

  for k = 1 : K
              
     w  = P(s(:,2),k);
     
E    = exp(egrid(k) + phip.*cond + phiu.*(~cond));               % compute tomorrow's dividends, given the optimal switching decision
rr   = funeval(cr,fspace,gridmake(x - kappap.*cond, k));         % tomorrow's multiplier, given tomorrow's state after pay fixed costs

L    = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*E;
K    = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E;
Y    = E.^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D    = Y - W.*L - (r+delta).*K;

     
     v  = v  + 1/(1+r)*w.*(funeval(c(:,3), fspace, gridmake(x, k)).*(~cond) + funeval(cp(:,3), fspacep, gridmake(x - kappap*cond, k)).*cond + D);       
     
  end