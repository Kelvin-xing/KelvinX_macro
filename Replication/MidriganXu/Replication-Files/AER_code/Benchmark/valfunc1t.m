  function v = valfunc1t(c,fspace,s,x)

  global beta
          
        v  = feval('menufun','ft',s,x) + beta*funeval(c(:,2),fspace,[x, s(:,2)]);
        