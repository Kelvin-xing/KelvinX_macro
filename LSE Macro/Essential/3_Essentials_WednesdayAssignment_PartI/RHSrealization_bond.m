% function output=RHSrealization_bond(epsi);
% Compute the actual RHS for realization shock
function output=RHSrealization_bond(epsi)

global T1 T po_k po_t dfactor nu rho cs lnts ks_S psi maturity iter psi_bond

t_next   = exp(rho*lnts(T1:T-1)+epsi); %note that t_next does not contain the first T1 observations
t_next_S = log(t_next);

xx       = makepoly([po_k po_t],[ks_S(T1:T-1) t_next_S]);

mu_next  = exp(xx*psi);
%note that mu_next does not contain the first T1 observations

%XXX 
%
%finish this function. Remember the thing to come up with at the end of the
%function (called "output") is the thing inside the conditional expectation
%for a realization of epsilon


if maturity < 2
    output     = 
else
    maturity_next = maturity-1;
    psi_bb  = psi_bond(:,iter,maturity_next,2);
    q_next  =                 ; %next period's bond price 
% it is possible to get the realization without specifying q-next. In fact 
% that may be more accurate since you don't use another approximation. But
% then you have to be careful with the sample used
    output  =                 ;
end


% **********************************************************************
% **********************************************************************

