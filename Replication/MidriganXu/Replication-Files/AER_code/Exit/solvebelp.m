function [bel, beljac] = solvebelp(c,fspace,s)

ns = size(s,1);

c  = [c(1:ns,:), c(ns+1:2*ns,:), c(2*ns+1:3*ns,:)];

[v1, v2, v3, ~, v12, v21, v33] = saveBelmaxp(c,fspace,s);

LHS = funeval(c,fspace,s);                                 % left-hand side of Bellman equation
B   = funbas(fspace,s);                                      % Miranda calls this the PHI matrix
                                                             
bel = zeros(3*length(s),1);

bel(1:ns)         =  LHS(:,1) - v1;                                   
bel(ns+1:2*ns)    =  LHS(:,2) - v2;   
bel(2*ns+1:3*ns)  =  LHS(:,3) - v3;                                 


v11 = zeros(ns,ns);
v22 = zeros(ns,ns); 
v13 = zeros(ns,ns);
v23 = zeros(ns,ns);
v31 = zeros(ns,ns);
v32 = zeros(ns,ns);


beljac = blkdiag(B,B, B) - [v11, v12, v13;...
                            v21, v22, v23;
                            v31, v32, v33];  
