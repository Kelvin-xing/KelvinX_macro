function [bel, beljac] = solvebelw(c,fspace,s)

[v, ~, vd] = saveBelmaxw(c, fspace, s);

LHS = funeval(c,fspace,s);                                 % left-hand side of Bellman equation
B   = funbas(fspace,s);                                      % Miranda calls this the PHI matrix
                                                             
bel    = LHS - v;
beljac = B   - vd; 