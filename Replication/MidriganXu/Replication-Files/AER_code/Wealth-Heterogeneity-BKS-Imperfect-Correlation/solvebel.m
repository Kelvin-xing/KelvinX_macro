function [bel, beljac] = solvebel(c,fspace,s)

ns = size(s,1);

c  = [c(1:ns,:), c(ns+1:2*ns,:), c(2*ns+1:3*ns,:), c(3*ns+1:4*ns,:)];

[v1, v2, v3, v4, ~, ~, v12, v21, v23, v34, v41, v43] = saveBelmax(c,fspace,s);

LHS = funeval(c,fspace,s);                                 % left-hand side of Bellman equation
B   = funbas(fspace,s);                                    % Miranda calls this the PHI matrix
                                                             
bel = zeros(4*length(s),1);

bel(1:ns)         =  LHS(:,1) - v1;                                   
bel(ns+1:2*ns)    =  LHS(:,2) - v2;   
bel(2*ns+1:3*ns)  =  LHS(:,3) - v3;   
bel(3*ns+1:4*ns)  =  LHS(:,4) - v4;   


v11 = zeros(ns,ns); v13 = zeros(ns,ns); v14 = zeros(ns,ns);
v22 = zeros(ns,ns); v24 = zeros(ns,ns); 
v31 = zeros(ns,ns); v32 = zeros(ns,ns); v33 = zeros(ns,ns); 
v42 = zeros(ns,ns); v44 = zeros(ns,ns); 

beljac = blkdiag(B, B, B, B) - [v11, v12, v13, v14;...
                                v21, v22, v23, v24;...
                                v31, v32, v33, v34;...
                                v41, v42, v43, v44];  