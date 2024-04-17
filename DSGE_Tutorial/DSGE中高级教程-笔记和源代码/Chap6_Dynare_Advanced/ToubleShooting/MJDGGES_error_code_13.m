%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'MJDGGES_error_code_13';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('MJDGGES_error_code_13.log');
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names_long = 'eps_a';
M_.exo_names = char(M_.exo_names, 'eps_y_star');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_y\_star');
M_.exo_names_long = char(M_.exo_names_long, 'eps_y_star');
M_.endo_names = 'y_star';
M_.endo_names_tex = 'y\_star';
M_.endo_names_long = 'y_star';
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'ygap');
M_.endo_names_tex = char(M_.endo_names_tex, 'ygap');
M_.endo_names_long = char(M_.endo_names_long, 'ygap');
M_.endo_names = char(M_.endo_names, 'rnat');
M_.endo_names_tex = char(M_.endo_names_tex, 'rnat');
M_.endo_names_long = char(M_.endo_names_long, 'rnat');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'ynat');
M_.endo_names_tex = char(M_.endo_names_tex, 'ynat');
M_.endo_names_long = char(M_.endo_names_long, 'ynat');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'pi_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_h');
M_.endo_names_long = char(M_.endo_names_long, 'pi_h');
M_.endo_names = char(M_.endo_names, 'pi_star');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_star');
M_.endo_names_long = char(M_.endo_names_long, 'pi_star');
M_.endo_names = char(M_.endo_names, 's');
M_.endo_names_tex = char(M_.endo_names_tex, 's');
M_.endo_names_long = char(M_.endo_names_long, 's');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'e');
M_.endo_names_tex = char(M_.endo_names_tex, 'e');
M_.endo_names_long = char(M_.endo_names_long, 'e');
M_.endo_names = char(M_.endo_names, 'p_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'p\_h');
M_.endo_names_long = char(M_.endo_names_long, 'p_h');
M_.endo_names = char(M_.endo_names, 'cpi_level');
M_.endo_names_tex = char(M_.endo_names_tex, 'cpi\_level');
M_.endo_names_long = char(M_.endo_names_long, 'cpi_level');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'nx');
M_.endo_names_tex = char(M_.endo_names_tex, 'nx');
M_.endo_names_long = char(M_.endo_names_long, 'nx');
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names_long = char(M_.param_names_long, 'sigma');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'rho_y_star');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_y\_star');
M_.param_names_long = char(M_.param_names_long, 'rho_y_star');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pi');
M_.param_names_long = char(M_.param_names_long, 'phi_pi');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_y');
M_.param_names_long = char(M_.param_names_long, 'phi_y');
M_.param_names = char(M_.param_names, 'kappa_a');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_a');
M_.param_names_long = char(M_.param_names_long, 'kappa_a');
M_.param_names = char(M_.param_names, 'omega');
M_.param_names_tex = char(M_.param_names_tex, 'omega');
M_.param_names_long = char(M_.param_names_long, 'omega');
M_.param_names = char(M_.param_names, 'sigma_a');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_a');
M_.param_names_long = char(M_.param_names_long, 'sigma_a');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, 'lambda');
M_.param_names_long = char(M_.param_names_long, 'lambda');
M_.param_names = char(M_.param_names, 'BigGamma_a');
M_.param_names_tex = char(M_.param_names_tex, 'BigGamma\_a');
M_.param_names_long = char(M_.param_names_long, 'BigGamma_a');
M_.param_names = char(M_.param_names, 'BigGamma_star');
M_.param_names_tex = char(M_.param_names_tex, 'BigGamma\_star');
M_.param_names_long = char(M_.param_names_long, 'BigGamma_star');
M_.param_names = char(M_.param_names, 'BigTheta');
M_.param_names_tex = char(M_.param_names_tex, 'BigTheta');
M_.param_names_long = char(M_.param_names_long, 'BigTheta');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 17;
M_.param_nbr = 18;
M_.orig_endo_nbr = 17;
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
erase_compiled_function('MJDGGES_error_code_13_static');
erase_compiled_function('MJDGGES_error_code_13_dynamic');
M_.lead_lag_incidence = [
 1 7 24;
 2 8 0;
 0 9 25;
 0 10 0;
 0 11 0;
 0 12 26;
 0 13 0;
 0 14 27;
 0 15 28;
 0 16 0;
 3 17 0;
 0 18 0;
 4 19 0;
 5 20 0;
 6 21 0;
 0 22 0;
 0 23 0;]';
M_.nstatic = 7;
M_.nfwrd   = 4;
M_.npred   = 5;
M_.nboth   = 1;
M_.nsfwrd   = 5;
M_.nspred   = 6;
M_.ndynamic   = 10;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(17, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(18, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 54;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.99;
beta = M_.params( 1 );
M_.params( 2 ) = 1;
sigma = M_.params( 2 );
M_.params( 3 ) = 0.4;
alpha = M_.params( 3 );
M_.params( 4 ) = 1;
eta = M_.params( 4 );
M_.params( 5 ) = 6;
epsilon = M_.params( 5 );
M_.params( 6 ) = 3;
phi = M_.params( 6 );
M_.params( 7 ) = 0.75;
theta = M_.params( 7 );
gamma =1; 
M_.params( 15 ) = (1-M_.params(1)*M_.params(7))*(1-M_.params(7))/M_.params(7);
lambda = M_.params( 15 );
M_.params( 13 ) = M_.params(2)*gamma+(1-M_.params(3))*(M_.params(2)*M_.params(4)-1);
omega = M_.params( 13 );
M_.params( 14 ) = M_.params(2)/(1+M_.params(3)*(M_.params(13)-1));
sigma_a = M_.params( 14 );
M_.params( 12 ) = M_.params(15)*(M_.params(6)+M_.params(14));
kappa_a = M_.params( 12 );
M_.params( 18 ) = M_.params(13)-1;
BigTheta = M_.params( 18 );
M_.params( 16 ) = (1+M_.params(6))/(M_.params(6)+M_.params(14));
BigGamma_a = M_.params( 16 );
M_.params( 17 ) = M_.params(14)*(-M_.params(3))*M_.params(18)/(M_.params(6)+M_.params(14));
BigGamma_star = M_.params( 17 );
M_.params( 8 ) = 0.66;
rho_a = M_.params( 8 );
M_.params( 9 ) = 0.86;
rho_y_star = M_.params( 9 );
M_.params( 10 ) = 1.5;
phi_pi = M_.params( 10 );
M_.params( 11 ) = 0;
phi_y = M_.params( 11 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
options_.irf = 16;
var_list_=[];
var_list_ = 'ygap';
var_list_ = char(var_list_, 'pi_h');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'pi');
var_list_ = char(var_list_, 'e');
var_list_ = char(var_list_, 'q');
var_list_ = char(var_list_, 'p_h');
var_list_ = char(var_list_, 'cpi_level');
info = stoch_simul(var_list_);
save('MJDGGES_error_code_13_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('MJDGGES_error_code_13_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('MJDGGES_error_code_13_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('MJDGGES_error_code_13_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('MJDGGES_error_code_13_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
