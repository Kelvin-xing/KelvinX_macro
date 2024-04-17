%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'cgg_level';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('cgg_level.log');
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names_long = 'eps_a';
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'f');
M_.endo_names_tex = char(M_.endo_names_tex, 'f');
M_.endo_names_long = char(M_.endo_names_long, 'f');
M_.param_names = 'gamma';
M_.param_names_tex = 'gamma';
M_.param_names_long = 'gamma';
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'Veps');
M_.param_names_tex = char(M_.param_names_tex, 'Veps');
M_.param_names_long = char(M_.param_names_long, 'Veps');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names_long = char(M_.param_names_long, 'cs');
M_.param_names = char(M_.param_names, 'ks');
M_.param_names_tex = char(M_.param_names_tex, 'ks');
M_.param_names_long = char(M_.param_names_long, 'ks');
M_.param_names = char(M_.param_names, 'as');
M_.param_names_tex = char(M_.param_names_tex, 'as');
M_.param_names_long = char(M_.param_names_long, 'as');
M_.param_names = char(M_.param_names, 'fs');
M_.param_names_tex = char(M_.param_names_tex, 'fs');
M_.param_names_long = char(M_.param_names_long, 'fs');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 10;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('cgg_level_static');
erase_compiled_function('cgg_level_dynamic');
M_.lead_lag_incidence = [
 0 3 7;
 1 4 0;
 2 5 8;
 0 6 0;]';
M_.nstatic = 1;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(10, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 13;
M_.NNZDerivatives(2) = 18;
M_.NNZDerivatives(3) = -1;
M_.params( 4 ) = 0.99;
beta = M_.params( 4 );
M_.params( 1 ) = 2;
gamma = M_.params( 1 );
M_.params( 2 ) = 0.36;
alpha = M_.params( 2 );
M_.params( 3 ) = 0.02;
delta = M_.params( 3 );
M_.params( 5 ) = 0.95;
rho = M_.params( 5 );
M_.params( 6 ) = 0.0001;
Veps = M_.params( 6 );
M_.params( 8 ) = (M_.params(2)*M_.params(4)/(1-M_.params(4)*(1-M_.params(3))))^(1/(1-M_.params(2)));
ks = M_.params( 8 );
M_.params( 9 ) = 0;
as = M_.params( 9 );
M_.params( 10 ) = M_.params(8)^M_.params(2)+(1-M_.params(3))*M_.params(8);
fs = M_.params( 10 );
M_.params( 7 ) = M_.params(10)-M_.params(8);
cs = M_.params( 7 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 3 ) = M_.params(9);
oo_.steady_state( 4 ) = log(M_.params(10));
oo_.steady_state( 2 ) = log(M_.params(8));
oo_.steady_state( 1 ) = log(M_.params(7));
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(6);
options_.nograph = 1;
options_.order = 2;
options_.qz_zero_threshold = 1e-15;
var_list_=[];
info = stoch_simul(var_list_);
save('cgg_level_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('cgg_level_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('cgg_level_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('cgg_level_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('cgg_level_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
