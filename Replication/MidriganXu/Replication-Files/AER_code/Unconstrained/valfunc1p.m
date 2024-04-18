  function v = valfunc1p(c,fspace,s,x)

  global beta
          
        v  = feval('menufun','fp',s,x) + beta*funeval(c(:,2),fspace,[x, s(:,2)]);
        