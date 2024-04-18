function [v1, v2, v3, v4, xm, xt, v12, v21, v23, v34, v41, v43] = saveBelmax(c,fspace,s)

global beta wnode znode rhom rhot kappa smin

ns = size(s,1);

solvem;

xm  = x; 
v1  = valfunc1m(c, fspace, s, xm);
v2  = valfunc2m(c, fspace, s);

solvet; 

xt  = x;
v3  = valfunc1t(c, fspace, s, xt);
v4  = valfunc2t(c, fspace, s);

if nargout > 6
  
v12 =  beta*funbas(fspace, [xm, s(:,2)]);
v34 =  beta*funbas(fspace, [xt, s(:,2)]);

v21  = 0;
v23  = 0;
v41  = 0; 
v43  = 0;  % if stay in trad., nothing changes 
        

vtt  = 0; 
vtm  = 0; 

K = size(wnode,1);

      for k = 1 : K
        
        w  = wnode(k);

        
        gs  = gridmake(s(:,1), znode(k));

        v43  = v43 + w*(rhom*funbas(fspace,s) + (1-rhom)*funbas(fspace,gs));
        
        vtt  = vtt + w*(rhom*funeval(c(:,3), fspace, s) + (1-rhom)*funeval(c(:,3), fspace, gs));  % if continue in trad., no changes: keep assets and efficiency
     
        v21 = v21 + w*rhom*funbas(fspace,s);
        
        g  = gridmake(s(:,1), znode(k));  

        v23 = v23 + w*(1-rhom)*funbas(fspace,g);
              
        g1  = [s(:,1) - kappa*exp(s(:,2)), s(:,2)];  
        g2  = gridmake(s(:,1) - kappa*exp(s(:,2)), znode(k));  
        
        vtm = vtm + w*(rhot*funeval(c(:,1), fspace, g1) + (1-rhot)*funeval(c(:,1), fspace, g2));
             
        v41 = v41 + w*(rhot*funbas(fspace,g1) + (1-rhot)*funbas(fspace,g2)); 
        
      end
     
     cond = vtm>=vtt & s(:,1) - kappa*exp(s(:,2)) > smin(1);    % can't move unless enough assets
     
     v41  = bsxfun(@times, v41, cond);
     v43  = bsxfun(@times, v43, ~cond);
      
end
