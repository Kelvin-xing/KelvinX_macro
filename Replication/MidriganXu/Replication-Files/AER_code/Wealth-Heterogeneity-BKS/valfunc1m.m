  function v = valfunc1m(c, fspace, s, x)

  global beta
          
        v  = feval('menufun','fm',s,x) + beta*funeval(c(:,2), fspace, [x, s(:,2)]);
        