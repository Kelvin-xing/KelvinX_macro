% function output=RHSrealization_R(epsi);
% Compute the actual RHS for realization shock
function output=RHSrealization_K(epsi)

global T1 T po_k po_t dfactor nu alpha rho delta cs lnts ks ks_S ts_mean ts_std psi

t_next   = exp(rho*lnts(T1:T-1)+epsi);
t_next_S = (t_next-ts_mean)/ts_std;

xx       = makepoly([po_k po_t],[ks_S(T1:T-1) t_next_S]);

%XXX
%
%finish this function. Remember the thing to come up with at the end of the
%function (called "output") is the thing inside the conditional expectation
%for a realization of epsilon


k_next  = ;
c_next  = ;
output     = ;


% **********************************************************************
% **********************************************************************

