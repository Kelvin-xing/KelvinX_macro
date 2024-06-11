//**************************************************************************
//                   LSE Macroeconomics Summer Program
//                   Part II: Heterogeneous Agents
//                   Instructor: Wouter J. Den Haan
//
//                   use of this program in any fee-based program requires
//                   explicit permission (wjdenhaan@gmail.com)
//**************************************************************************
//
%==========================================================================
%                      Simple Aiyagari type model
%==========================================================================

% This file solves the Aiyagari type model for a given value of the 
% aggregate interest rate (implying a value for aggregate capital demand 
% and the aggregate wage rate). Instead of a non-borrowing constraint we 
% employ a penalty function on capital accumulation. 

%--------------------------------------------------------------------------
% Declarations
%--------------------------------------------------------------------------

var c, k, z;
varexo e;
parameters betta, alpha, delta, zeta0, zeta1, zeta2, r, w, rho, 
       sigshock, nu, c_ss, k_init;

%--------------------------------------------------------------------------
% 2. Parameter values
%--------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below load and set all the necessary parameter values
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

load ;
set_param_value


kd = (r/alpha)^(1/(alpha-1));
w  = (1-alpha)*kd^alpha;
c_ss = ((kd)*(r-delta) + w);

%--------------------------------------------------------------------------
% 3. Model equations
%--------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below fill in the model block
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

model;


end;

%--------------------------------------------------------------------------
% 4. Steady states
%--------------------------------------------------------------------------

initval;

k = k_init;
c = c_ss;
z = 1;

end;

%--------------------------------------------------------------------------
% 5. Shocks and solution
%--------------------------------------------------------------------------

shocks;
var e; stderr sigshock;
end;

steady;
stoch_simul(order=2,nocorr,nomoments,IRF=0) k;