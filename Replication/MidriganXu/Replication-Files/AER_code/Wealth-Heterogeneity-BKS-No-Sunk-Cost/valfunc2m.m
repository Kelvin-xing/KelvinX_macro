  function v = valfunc2m(c,fspace,s)

  global wnode znode rhom
        
K = size(wnode,1);

v = 0;

      for k = 1 : K
        
        w  = wnode(k);

        g  = gridmake(s(:,1), znode(k));  
        
        v = v + w.*(rhom*funeval(c(:,1), fspace, s) + (1-rhom)*funeval(c(:,3), fspace, g));
                      
      end

   