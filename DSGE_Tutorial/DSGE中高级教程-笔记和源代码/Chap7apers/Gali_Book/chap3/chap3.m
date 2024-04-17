%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'chap3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('chap3.log');
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names_long = 'eps_a';
M_.exo_names = char(M_.exo_names, 'eps_nu');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_nu');
M_.exo_names_long = char(M_.exo_names_long, 'eps_nu');
M_.endo_names = 'pi';
M_.endo_names_tex = 'pi';
M_.endo_names_long = 'pi';
M_.endo_names = char(M_.endo_names, 'y_gap');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_gap');
M_.endo_names_long = char(M_.endo_names_long, 'y_gap');
M_.endo_names = char(M_.endo_names, 'y_nat');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_nat');
M_.endo_names_long = char(M_.endo_names_long, 'y_nat');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'r_nat');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_nat');
M_.endo_names_long = char(M_.endo_names_long, 'r_nat');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'm_growth_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'm\_growth\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'm_growth_ann');
M_.endo_names = char(M_.endo_names, 'nu');
M_.endo_names_tex = char(M_.endo_names_tex, 'nu');
M_.endo_names_long = char(M_.endo_names_long, 'nu');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'r_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'r_ann');
M_.endo_names = char(M_.endo_names, 'R_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'R\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'R_ann');
M_.endo_names = char(M_.endo_names, 'r_nat_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_nat\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'r_nat_ann');
M_.endo_names = char(M_.endo_names, 'pi_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'pi_ann');
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names_long = char(M_.param_names_long, 'sigma');
M_.param_names = char(M_.param_names, 'psi_n_ya');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_n\_ya');
M_.param_names_long = char(M_.param_names_long, 'psi_n_ya');
M_.param_names = char(M_.param_names, 'rho_nu');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_nu');
M_.param_names_long = char(M_.param_names_long, 'rho_nu');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pi');
M_.param_names_long = char(M_.param_names_long, 'phi_pi');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_y');
M_.param_names_long = char(M_.param_names_long, 'phi_y');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 15;
M_.param_nbr = 12;
M_.orig_endo_nbr = 15;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('chap3_static');
erase_compiled_function('chap3_dynamic');
M_.lead_lag_incidence = [
 0 5 20;
 0 6 21;
 0 7 0;
 1 8 0;
 0 9 0;
 0 10 0;
 2 11 0;
 0 12 0;
 0 13 0;
 3 14 0;
 4 15 22;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 0;]';
M_.nstatic = 9;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(15, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 46;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 2 ) = 1;
sigma = M_.params( 2 );
M_.params( 6 ) = 1;
phi = M_.params( 6 );
M_.params( 7 ) = 1.5;
phi_pi = M_.params( 7 );
M_.params( 8 ) = 0.125;
phi_y = M_.params( 8 );
theta=2/3;      
M_.params( 4 ) = 0.5;
rho_nu = M_.params( 4 );
M_.params( 5 ) = 0.9;
rho_a = M_.params( 5 );
M_.params( 1 ) = 0.99;
beta = M_.params( 1 );
M_.params( 12 ) = 4;
eta = M_.params( 12 );
M_.params( 10 ) = 0.3333333333333333;
alpha = M_.params( 10 );
M_.params( 11 ) = 6;
epsilon = M_.params( 11 );
Omega=(1-alpha)/(1-alpha+alpha*epsilon);  
M_.params( 3 ) = (1+M_.params(6))/(M_.params(10)+M_.params(6)+M_.params(2)*(1-M_.params(10)));
psi_n_ya = M_.params( 3 );
lambda=(1-theta)*(1-beta*theta)/theta*Omega; 
M_.params( 9 ) = lambda*(M_.params(2)+(M_.params(6)+M_.params(10))/(1-M_.params(10)));
kappa = M_.params( 9 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(2, 2) = 0.0625;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 15;
options_.order = 1;
var_list_=[];
var_list_ = 'y_gap';
var_list_ = char(var_list_, 'pi_ann');
var_list_ = char(var_list_, 'R_ann');
var_list_ = char(var_list_, 'r_ann');
var_list_ = char(var_list_, 'm_growth_ann');
var_list_ = char(var_list_, 'nu');
info = stoch_simul(var_list_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 0;
options_.irf = 15;
options_.order = 1;
var_list_=[];
var_list_ = 'y_gap';
var_list_ = char(var_list_, 'pi_ann');
var_list_ = char(var_list_, 'y');
var_list_ = char(var_list_, 'n');
var_list_ = char(var_list_, 'R_ann');
var_list_ = char(var_list_, 'r_ann');
var_list_ = char(var_list_, 'm_growth_ann');
var_list_ = char(var_list_, 'a');
info = stoch_simul(var_list_);
save('chap3_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('chap3_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('chap3_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('chap3_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('chap3_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
