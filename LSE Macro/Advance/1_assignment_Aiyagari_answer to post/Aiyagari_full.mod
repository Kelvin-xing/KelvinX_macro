%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
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


load parametervalues;
set_param_value('sigshock'   ,par.sigshock);
set_param_value('r'          ,r);
set_param_value('betta'      ,par.beta);
set_param_value('alpha'      ,par.alpha);
set_param_value('delta'      ,par.delta);
set_param_value('rho'        ,par.rho);
set_param_value('nu'         ,par.nu);

set_param_value('k_init'     ,par.kini);
set_param_value('zeta0'      ,par.zeta0);
set_param_value('zeta1'      ,par.zeta1);
set_param_value('zeta2'      ,par.zeta2);


kd = (r/alpha)^(1/(alpha-1));
w  = (1-alpha)*kd^alpha;
c_ss = ((kd)*(r-delta) + w);

%--------------------------------------------------------------------------
% 3. Model equations
%--------------------------------------------------------------------------

model;

c^(-nu) = -zeta2+zeta1*exp(-zeta0*k) + betta*c(+1)^(-nu)*(r + 1 - delta);
c + k   = r*k(-1) + w*z + (1-delta)*k(-1);
z = 1-rho + rho*z(-1) + e;
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
stoch_simul(order=2,nocorr,nomoments,IRF=0);