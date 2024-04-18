  function v = valfunc2u(c,fspace,s)

  global P
        
K = size(P,2);

v = 0;

      for k = 1 : K
        
        w  = P(s(:,2),k);

        g  = gridmake(s(:,1), k);           
        v = v + w.*funeval(c(:,1), fspace, g);
                      
      end

   