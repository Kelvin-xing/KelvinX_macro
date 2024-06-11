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

var z_1      , z_2      , z_i      , 
    z_1_exp_A, z_2_exp_A, z_i_exp_A,
    n_1_R, n_2_R, y_R, 
    n_1_A, n_2_A, y_A, 
    wage, wage_exp_R, wage_exp_A, 
    N,       N_exp_R,    N_exp_A;

varexo e_1,e_2,e_i;

parameters rho_z  , rho_i  , sig_z, sig_i, alpha, 
           rho_z_A, rho_i_A, eta,
           coef_N_z1  , coef_N_z2  , coef_N_Nlag  ,
           coef_N_z1_A, coef_N_z2_A, coef_N_Nlag_A,
           wcoef, wcoef_A, betta;


load parametervalues;
set_param_value('alpha'    ,alpha)
set_param_value('rho_z'    ,rho_z)
set_param_value('rho_z_A'  ,rho_z_A)
set_param_value('rho_i'    ,rho_i)
set_param_value('rho_i_A'  ,rho_i_A)
set_param_value('sig_z'    ,sig_z)
set_param_value('sig_i'    ,sig_i)
set_param_value('eta'      ,eta)
set_param_value('wcoef'    ,wcoef)
set_param_value('wcoef_A'  ,wcoef_A)
set_param_value('betta'  ,betta)


load aggregatelaw_N;
set_param_value('coef_N_z1'   ,  coef_N_z1)
set_param_value('coef_N_z2'     ,coef_N_z2)
set_param_value('coef_N_Nlag'   ,coef_N_Nlag)
set_param_value('coef_N_z1_A'   ,coef_N_z1_A)
set_param_value('coef_N_z2_A'   ,coef_N_z2_A)
set_param_value('coef_N_Nlag_A' ,coef_N_Nlag_A)

model;

%true laws of motion
z_1 = 1-rho_z   + rho_z  * z_1(-1) + e_1;
z_2 = 1-rho_z   + rho_z  * z_2(-1) + e_2;
z_i = 1-rho_i   + rho_i  * z_i(-1) + e_i;

%expectations according to perceived laws of motion by type A firms
z_1_exp_A = 1-rho_z_A  + rho_z_A * z_1;
z_2_exp_A = 1-rho_z_A  + rho_z_A * z_2;
z_i_exp_A = 1-rho_i_A  + rho_i_A * z_i;

% production levels
y_A = z_i*z_1*n_1_A^alpha + z_i*z_2*n_2_A^alpha;
y_R = z_i*z_1*n_1_R^alpha + z_i*z_2*n_2_R^alpha;

% FOC conditions
betta*alpha*n_1_R^(alpha-1)*z_i(+1)*z_1(+1) = betta*wage_exp_R + eta*(n_1_R    +n_2_R    - n_1_R(-1)-n_2_R(-1)) 
                                            - betta*eta*(n_1_R(+1)+n_2_R(+1)- n_1_R    -n_2_R    );
betta*alpha*n_2_R^(alpha-1)*z_i(+1)*z_2(+1) = betta*wage_exp_R + eta*(n_1_R+n_2_R- n_1_R(-1)-n_2_R(-1))
                                            - betta*eta*(n_1_R(+1)+n_2_R(+1)- n_1_R    -n_2_R    );

betta*alpha*n_1_A^(alpha-1)*z_1_exp_A*z_i_exp_A = betta*wage_exp_A + eta*(n_1_A    +n_2_A    - n_1_A(-1)-n_2_A(-1))
                                                - betta*eta*(n_1_A(+1)+n_2_A(+1)- n_1_A    -n_2_A    );
betta*alpha*n_2_A^(alpha-1)*z_2_exp_A*z_i_exp_A = betta*wage_exp_A + eta*(n_1_A    +n_2_A    - n_1_A(-1)-n_2_A(-1))
                                                - betta*eta*(n_1_A(+1)+n_2_A(+1)- n_1_A    -n_2_A    );


% true law of motion
N  = coef_N_z1*(z_1-1) + coef_N_z2*(z_2-1) + coef_N_Nlag*N(-1);

% perceived laws of motion
N_exp_R = coef_N_z1  *(z_1-1) + coef_N_z2  *(z_2-1) + coef_N_Nlag  *N(-1);
N_exp_A = coef_N_z1_A*(z_1-1) + coef_N_z2_A*(z_2-1) + coef_N_Nlag_A*N(-1);

% true wage rate
wage = alpha*0.5^(alpha-1) + wcoef*N(-1);

% expected wage rates
wage_exp_R = alpha*0.5^(alpha-1) + wcoef  *N_exp_R;
wage_exp_A = alpha*0.5^(alpha-1) + wcoef_A*N_exp_A;

end;

initval;
n_1_R     =		 0.5;
n_2_R     =      1-n_1_R;
n_1_A     =      0.5;
n_2_A     =      1-n_1_A;
y_R       =      2*0.5^alpha;
y_A       =      y_R;
z_1       =      1;
z_2       =      1;
z_1_exp_A =      1;
z_2_exp_A =      1;
z_i       =      1;
z_i_exp_A =      1;
N         =      0;
wage      = alpha*0.5^(alpha-1);
end;

shocks;
var e_1; stderr sig_z;
var e_2; stderr sig_z;
var e_i; stderr sig_i;
end;

steady;
stoch_simul(order=1,nocorr,nomoments,IRF=0) n_1_R n_2_R n_1_A n_2_A N;
