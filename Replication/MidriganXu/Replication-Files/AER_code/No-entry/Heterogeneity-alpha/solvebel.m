function [bel, beljac] = solvebel(c,fspace,s)

ns = size(s,1);

c  = [c(1:ns,:), c(ns+1:2*ns,:)];

[v1, v2, ~, v12, v21] = saveBelmax(c,fspace,s);

LHS = funeval(c,fspace,s);                                 % left-hand side of Bellman equation
B   = funbas(fspace,s);                                      % Miranda calls this the PHI matrix
                                                             
bel = zeros(2*length(s),1);

bel(1:ns)         =  LHS(:,1) - v1;                                   
bel(ns+1:2*ns)    =  LHS(:,2) - v2;   


v11 = zeros(ns,ns);
v22 = zeros(ns,ns); 

beljac = blkdiag(B,B) - [v11, v12;...
                         v21, v22];  
