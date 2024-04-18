  function v = valfunc1u(c,fspace,s,x)

  global beta
          
        v  = feval('menufun','fu',s,x) + beta*funeval(c(:,2),fspace,[x, s(:,2)]);
        