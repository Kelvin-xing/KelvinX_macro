  function v = valfunc2(c,fspace,s)

  global P
  
        v  = 0;
        ns = size(s,1);
        
K = size(P,2);
     
  for k = 1 : K
              
     w  = P(s(:,2),k);
     
     v  = v  + w.*funeval(c(:,1), fspace, [s(:,1), k*ones(ns,1), s(:,2)]);       
     
  end

  
