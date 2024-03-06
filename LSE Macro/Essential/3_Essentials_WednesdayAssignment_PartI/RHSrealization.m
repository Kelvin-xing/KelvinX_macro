% function output=RHSrealization(epsi);
% Compute the RHS for realization shock
function output=RHSrealization(epsi)

global T1 T po_k po_t dfactor nu alpha rho delta lnts ks ks_S psi

% calculate next period's TFP for next period's realization of shock epsi
% lnts(t) is the (log of) TFP in period t

t_next   = exp(rho*lnts(T1:T-1)+epsi); %note that t_next does not contain the first T1 observations

t_next_S = log(t_next);

% calculate polynomial basis functions using next period's state variables
% ks_S(t) is the (scaled) level of capital chosen in period T1 and thus
% state variable in period t+1s

polynomial  = makepoly([po_k po_t],[ks_S(T1:T-1) t_next_S]);
%note that you don't have exclude any observations in t_next_S

%XXX
%
%finish this function. Remember the thing to come up with at the end of the
%function (called "output") is the thing inside the conditional expectation
%for a realization of epsilon

output = 


% **********************************************************************
% **********************************************************************

