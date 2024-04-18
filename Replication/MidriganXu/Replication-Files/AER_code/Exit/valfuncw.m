  function v = valfuncw(c,fspace,s,x)

  global beta Pw
          
  K = size(Pw, 2);
     
  v  = feval('menufunw','fw',s,x);
  
  for k = 1 : K
              
     w  = Pw(s(:,2),k);

     v  = v  + beta*w.*funeval(c, fspace, gridmake(x, k));       
     
  end

  