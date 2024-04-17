%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'bgg_rbc';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('bgg_rbc.log');
M_.exo_names = 'sigma_e';
M_.exo_names_tex = 'sigma\_e';
M_.exo_names_long = 'sigma_e';
M_.endo_names = 'kbar';
M_.endo_names_tex = 'kbar';
M_.endo_names_long = 'kbar';
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'omegabar');
M_.endo_names_tex = char(M_.endo_names_tex, 'omegabar');
M_.endo_names_long = char(M_.endo_names_long, 'omegabar');
M_.endo_names = char(M_.endo_names, 'Rk');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rk');
M_.endo_names_long = char(M_.endo_names_long, 'Rk');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'sigma');
M_.endo_names_tex = char(M_.endo_names_tex, 'sigma');
M_.endo_names_long = char(M_.endo_names_long, 'sigma');
M_.endo_names = char(M_.endo_names, 'spread');
M_.endo_names_tex = char(M_.endo_names_tex, 'spread');
M_.endo_names_long = char(M_.endo_names_long, 'spread');
M_.endo_names = char(M_.endo_names, 'credit');
M_.endo_names_tex = char(M_.endo_names_tex, 'credit');
M_.endo_names_long = char(M_.endo_names_long, 'credit');
M_.endo_names = char(M_.endo_names, 'bankrupt');
M_.endo_names_tex = char(M_.endo_names_tex, 'bankrupt');
M_.endo_names_long = char(M_.endo_names_long, 'bankrupt');
M_.endo_names = char(M_.endo_names, 'GDP');
M_.endo_names_tex = char(M_.endo_names_tex, 'GDP');
M_.endo_names_long = char(M_.endo_names_long, 'GDP');
M_.endo_names = char(M_.endo_names, 'wedge');
M_.endo_names_tex = char(M_.endo_names_tex, 'wedge');
M_.endo_names_long = char(M_.endo_names_long, 'wedge');
M_.param_names = 'sigma_ss';
M_.param_names_tex = 'sigma\_ss';
M_.param_names_long = 'sigma_ss';
M_.param_names = char(M_.param_names, 'mu');
M_.param_names_tex = char(M_.param_names_tex, 'mu');
M_.param_names_long = char(M_.param_names_long, 'mu');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'rhosigma');
M_.param_names_tex = char(M_.param_names_tex, 'rhosigma');
M_.param_names_long = char(M_.param_names_long, 'rhosigma');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, 'gam');
M_.param_names_long = char(M_.param_names_long, 'gam');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 13;
M_.param_nbr = 8;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('bgg_rbc_static');
erase_compiled_function('bgg_rbc_dynamic');
M_.lead_lag_incidence = [
 1 6 0;
 0 7 0;
 0 8 19;
 2 9 0;
 0 10 20;
 0 11 21;
 3 12 0;
 4 13 0;
 0 14 0;
 5 15 0;
 0 16 0;
 0 17 0;
 0 18 0;]';
M_.nstatic = 5;
M_.nfwrd   = 3;
M_.npred   = 5;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 5;
M_.ndynamic   = 8;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 50;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.26;
sigma_ss = M_.params( 1 );
M_.params( 2 ) = 0.21;
mu = M_.params( 2 );
M_.params( 3 ) = 0.97;
gamma = M_.params( 3 );
M_.params( 4 ) = 0.3333333333333333;
alpha = M_.params( 4 );
M_.params( 5 ) = 0.02;
delta = M_.params( 5 );
M_.params( 6 ) = 0.9926375361451395;
beta = M_.params( 6 );
M_.params( 7 ) = 0.97;
rhosigma = M_.params( 7 );
M_.params( 8 ) = 1;
gam = M_.params( 8 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (.11)^2;
steady;
options_.irf = 20;
options_.order = 1;
options_.periods = 1000;
var_list_=[];
var_list_ = 'kbar';
var_list_ = char(var_list_, 'i');
var_list_ = char(var_list_, 'c');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'omegabar');
var_list_ = char(var_list_, 'Rk');
var_list_ = char(var_list_, 'n');
var_list_ = char(var_list_, 'sigma');
var_list_ = char(var_list_, 'spread');
var_list_ = char(var_list_, 'credit');
var_list_ = char(var_list_, 'bankrupt');
var_list_ = char(var_list_, 'GDP');
var_list_ = char(var_list_, 'wedge');
info = stoch_simul(var_list_);
save('bgg_rbc_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('bgg_rbc_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('bgg_rbc_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('bgg_rbc_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('bgg_rbc_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
