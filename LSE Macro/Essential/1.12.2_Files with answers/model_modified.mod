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

var c a y p1 pn;                                                //Added: p1 and pn
varexo e;
parameters beta gamma mu sigma r eta0 eta1 eta2 a0 a1;          //Added: a0 and a1 (1st-order perturbation solution,
                                                                //already solved for by preliminary run of Dynare)

//  Parameters
//------------

//Includes parameters.mod which is created in main program.
@#include "parameters.mod"
a0 = ...;                                                       //Added: value for a0, but could also be included in "parameters.mod"
a1 = ...;                                                       //Added: value for a1, but could also be included in "parameters.mod"

//  Model
//-------

model;
//Core model block is unchanged
c^-gamma/(1+r) = beta*c(+1)^-gamma + eta1*exp(-eta0*a)+eta2;
c + a/(1+r) = a(-1) + y;
y = mu + e;
//Just include two extra equations for p1 and pn
p1 = a0 + a1 * (a(-1)-a0+e);					//Added: equation for p1 (1st-order perturbation solution)
a  = exp(-(a(-1)-a0+e)^2) * pn + (1-exp(-(a(-1)-a0+e)^2)) * p1; //Added: equation for pn (weighted average as in assignment)
end;

//Note: you can also include the correction term for p1.
//Then write: p1 = a0 + correction term + a1 * (a(-1)-a0+e).

//  Steady state
//--------------

initval;
//Includes steadystate.mod which is created in main program.
@#include "steadystate.mod"
p1 = a;                                                         //Added: initval for p1, but could also be included in "steadystate.mod"
pn = a;                                                         //Added: initval for pn, but could also be included in "steadystate.mod"
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