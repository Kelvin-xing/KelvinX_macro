//               Amsterdam Macroeconomics Summer School 2010
//                      Part II: Heterogeneous Agents
//                         University of Amsterdam
//
//                           Thursday Assignment
//            Badly behaved higher-order pertubation solutions &
//                        why pruning is a bad idea
//-------------------------------------------------------------------------

//  Declarations
//--------------

var c a y;
varexo e;
parameters beta gamma mu sigma r eta0 eta1 eta2;

//  Parameters
//------------

//Includes parameters.mod which is created in main program.
@#include "parameters.mod"

//  Model
//-------

model;
c^-gamma/(1+r) = beta*c(+1)^-gamma + eta1*exp(-eta0*a)+eta2;
c + a/(1+r) = a(-1) + y;
y = mu + e;
end;

//  Steady state
//--------------

initval;
//Includes steadystate.mod which is created in main program.
@#include "steadystate.mod"
end;

steady;

//  Shocks
//--------

shocks;
var e = sigma^2;
end;

//  Solution
//----------

stoch_simul(order=2,nocorr,nomoments,IRF=0);