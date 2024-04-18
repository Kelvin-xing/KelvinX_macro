  function v = valfunc2p(c,fspace,s)

  global P
  
        v  = 0;
        
K = size(P,2);
     
  for k = 1 : K
              
     w  = P(s(:,2),k);
     
     v  = v  + w.*funeval(c(:,1), fspace, gridmake(s(:,1), k));       
     
  end

  
