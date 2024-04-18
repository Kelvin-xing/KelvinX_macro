  function [v, vtu, vtt] = valfunc2t(c,fspace,s)

  global P fspaceu cu sminu kappau
 
     vtu  = 0;
     vtt  = 0;
        
K = size(P,2);

      for k = 1 : K
        
        w  = P(s(:,2),k);

        gt  = gridmake(s(:,1), k);           
        vtt = vtt + w.*funeval(c(:,1), fspace, gt);

        gu  = gridmake(s(:,1) - kappau, k);           
        vtu = vtu + w.*funeval(cu(:,1), fspaceu, gu);
                      
      end

     cond = vtu>=vtt & gu(:,1) > sminu(1);    % can't move unless sufficiently large stock of assets -- Amin
     
     v = vtu.*cond + vtt.*(~cond);
     
vtu(gu(:,1) < sminu(1)) = vtt(gu(:,1) < sminu(1)) - 5;