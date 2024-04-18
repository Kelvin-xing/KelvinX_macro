  function v = valfunc1(c,fspace,s,x)

  global beta
          
        v  = feval('menufun','f',s,x) + beta*funeval(c(:,2),fspace,[x, s(:,2)]);
        