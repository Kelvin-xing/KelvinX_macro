  function v = valfunc1t(c,fspace,s,x)

  global beta
          
        v  = feval('menufun','ft',s,x) + beta*funeval(c(:,4), fspace, [x, s(:,2)]);
        