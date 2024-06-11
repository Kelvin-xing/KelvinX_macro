%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
// heterobeliefs.mod
//

var z_1      , z_2      
    n_1_R, n_2_R, y_R, 
    wage, N;

varexo e_1,e_2;

%parameters rho_z  , rho_i  , sig_z, sig_i, alpha, 
%           rho_z_A, rho_i_A, eta,
%           coef_N_z1, coef_N_z2, coef_N_Nlag,
%           wcoef,betta;

parameters rho_z  , rho_i  , sig_z, sig_i, alpha, 
           eta,
           coef_N_z1, coef_N_z2, coef_N_Nlag,
           wcoef,betta;


load parametervalues;
set_param_value('alpha'    ,alpha)
set_param_value('rho_z'    ,rho_z)
set_param_value('rho_i'    ,rho_i)
set_param_value('sig_z'    ,sig_z)
set_param_value('sig_i'    ,sig_i)
set_param_value('eta'      ,eta)
set_param_value('wcoef'    ,wcoef)
set_param_value('betta'  ,betta)

model;
z_1 = 1-rho_z   + rho_z  * z_1(-1) + e_1;
z_2 = 1-rho_z   + rho_z  * z_2(-1) + e_2;

y_R = z_1*n_1_R^alpha + z_2*n_2_R^alpha;

betta*alpha*n_1_R^(alpha-1)*z_1(+1) = betta*wage(+1) + eta*(n_1_R    +n_2_R    - n_1_R(-1)-n_2_R(-1)) 
                                    - betta*eta*(n_1_R(+1)+n_2_R(+1)- n_1_R    -n_2_R    );
betta*alpha*n_2_R^(alpha-1)*z_2(+1) = betta*wage(+1) + eta*(n_1_R+n_2_R- n_1_R(-1)-n_2_R(-1))
                                    - betta*eta*(n_1_R(+1)+n_2_R(+1)- n_1_R    -n_2_R    );


wage = alpha*0.5^(alpha-1) + wcoef*(n_1_R(-1)+n_2_R(-1)-1);
N  = n_1_R + n_2_R - 1;
end;

initval;
n_1_R     =		 0.5;
n_2_R     =      1-n_1_R;
y_R       =      2*0.5^alpha;
z_1       =      1;
z_2       =      1;
N         =      0;
wage      = alpha*0.5^(alpha-1);
end;

shocks;
var e_1; stderr sig_z;
var e_2; stderr sig_z;
end;

steady;
//stoch_simul(order=1,solve_algo=3,nocorr,nomoments,replic=50,IRF=20) n_1_R n_1_A;
stoch_simul(order=1,nocorr,nomoments,IRF=20) n_1_R n_2_R N;
