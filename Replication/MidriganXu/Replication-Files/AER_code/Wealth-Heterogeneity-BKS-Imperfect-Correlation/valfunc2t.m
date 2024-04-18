  function [v, vtm, vtt] = valfunc2t(c, fspace, s)

  global wnode znode smin kappa rhot rhom

     vtt  = 0;
     vtm  = 0;
        
K = size(wnode, 1);

      for k = 1 : K
        
        w  = wnode(k);

        gs  = gridmake(s(:,1), znode(k));
        
        vtt  = vtt + w*(rhom*funeval(c(:,3), fspace, s) + (1-rhom)*funeval(c(:,3), fspace, gs));  % if continue in trad., no changes: keep assets and efficiency

        g1  = [s(:,1) - kappa*exp(s(:,2)), s(:,2)];  
        g2  = gridmake(s(:,1) - kappa*exp(s(:,2)), znode(k));  
        
        vtm = vtm + w*(rhot*funeval(c(:,1), fspace, g1) + (1-rhot)*funeval(c(:,1), fspace, g2));
                      
      end

     cond = vtm>=vtt & s(:,1) - kappa*exp(s(:,2)) > smin(1);    % can't move unless enough assets
     
     v = vtm.*cond + vtt.*(~cond);
    
     
vtm(s(:,1) - kappa*exp(s(:,2)) < smin(1)) = vtt(s(:,1) - kappa*exp(s(:,2)) < smin(1)) - 100;