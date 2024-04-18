a = sminu(1)*ones(size(si,1),1); 
b = smaxu(1)*ones(size(si,1),1);   % do not extrapolate

fmin = a - si(:,1) + kappau - theta*xai*funeval(cu(:,3), fspaceu, [a, si(:,2)]);
fmax = b - si(:,1) + kappau - theta*xai*funeval(cu(:,3), fspaceu, [b, si(:,2)]);

if any(fmin>0)
    
    b(fmin>0) = a(fmin>0);   % stop searching, can't enter with this level of savings. 
    
end

sb = sign(fmax);

% Initializations  

tol = 1e-11;

 dx = 0.5*(b - a);
 x  = a + dx;
 dx = sb.*dx;

% Iteration loop

  while any(abs(dx)>tol)

      dx = 0.5*dx;
      A = x;
      
sx = sign(A  - si(:,1) + kappau - theta*xai*funeval(cu(:,3), fspaceu, [A, si(:,2)]));   % this must be > 0
      
       x = x - sx.*dx;
     
  end

  
  
inject = theta*xai*funeval(cu(:,3), fspaceu, [A, si(:,2)]).*(fmin<=0);

