%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'indeterminacry_error';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('indeterminacry_error.log');
M_.exo_names = 'eps_A';
M_.exo_names_tex = 'eps\_A';
M_.exo_names_long = 'eps_A';
M_.exo_names = char(M_.exo_names, 'eps_m');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_m');
M_.exo_names_long = char(M_.exo_names_long, 'eps_m');
M_.endo_names = 'C';
M_.endo_names_tex = 'C';
M_.endo_names_long = 'C';
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.endo_names = char(M_.endo_names, 'N');
M_.endo_names_tex = char(M_.endo_names_tex, 'N');
M_.endo_names_long = char(M_.endo_names_long, 'N');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'm_growth_ann');
M_.endo_names_tex = char(M_.endo_names_tex, 'm\_growth\_ann');
M_.endo_names_long = char(M_.endo_names_long, 'm_growth_ann');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names_long = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names_long = char(M_.param_names_long, 'sigma');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pi');
M_.param_names_long = char(M_.param_names_long, 'phi_pi');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'Cs');
M_.param_names_tex = char(M_.param_names_tex, 'Cs');
M_.param_names_long = char(M_.param_names_long, 'Cs');
M_.param_names = char(M_.param_names, 'ws');
M_.param_names_tex = char(M_.param_names_tex, 'ws');
M_.param_names_long = char(M_.param_names_long, 'ws');
M_.param_names = char(M_.param_names, 'pis');
M_.param_names_tex = char(M_.param_names_tex, 'pis');
M_.param_names_long = char(M_.param_names_long, 'pis');
M_.param_names = char(M_.param_names, 'As');
M_.param_names_tex = char(M_.param_names_tex, 'As');
M_.param_names_long = char(M_.param_names_long, 'As');
M_.param_names = char(M_.param_names, 'Ns');
M_.param_names_tex = char(M_.param_names_tex, 'Ns');
M_.param_names_long = char(M_.param_names_long, 'Ns');
M_.param_names = char(M_.param_names, 'Rs');
M_.param_names_tex = char(M_.param_names_tex, 'Rs');
M_.param_names_long = char(M_.param_names_long, 'Rs');
M_.param_names = char(M_.param_names, 'rs');
M_.param_names_tex = char(M_.param_names_tex, 'rs');
M_.param_names_long = char(M_.param_names_long, 'rs');
M_.param_names = char(M_.param_names, 'Ys');
M_.param_names_tex = char(M_.param_names_tex, 'Ys');
M_.param_names_long = char(M_.param_names_long, 'Ys');
M_.param_names = char(M_.param_names, 'm_growth_anns');
M_.param_names_tex = char(M_.param_names_tex, 'm\_growth\_anns');
M_.param_names_long = char(M_.param_names_long, 'm_growth_anns');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 9;
M_.param_nbr = 16;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('indeterminacry_error_static');
erase_compiled_function('indeterminacry_error_dynamic');
M_.lead_lag_incidence = [
 0 4 13;
 0 5 0;
 0 6 14;
 1 7 0;
 0 8 0;
 2 9 0;
 0 10 0;
 3 11 0;
 0 12 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(16, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 30;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.33;
alpha = M_.params( 1 );
M_.params( 2 ) = 0.99;
beta = M_.params( 2 );
M_.params( 3 ) = 0.9;
rho = M_.params( 3 );
M_.params( 4 ) = 1;
sigma = M_.params( 4 );
M_.params( 5 ) = 1;
phi = M_.params( 5 );
M_.params( 6 ) = 1.5;
phi_pi = M_.params( 6 );
M_.params( 7 ) = 4;
eta = M_.params( 7 );
M_.params( 11 ) = 1;
As = M_.params( 11 );
M_.params( 13 ) = 1/M_.params(2);
Rs = M_.params( 13 );
M_.params( 10 ) = 1;
pis = M_.params( 10 );
M_.params( 14 ) = M_.params(13);
rs = M_.params( 14 );
M_.params( 12 ) = (1-M_.params(1))^(1/(M_.params(4)+M_.params(1)*(1-M_.params(4))+M_.params(5)));
Ns = M_.params( 12 );
M_.params( 8 ) = M_.params(11)*M_.params(12)^(1-M_.params(1));
Cs = M_.params( 8 );
M_.params( 9 ) = (1-M_.params(1))*M_.params(11)*M_.params(12)^(-M_.params(1));
ws = M_.params( 9 );
M_.params( 15 ) = M_.params(8);
Ys = M_.params( 15 );
M_.params( 16 ) = 0;
m_growth_anns = M_.params( 16 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = log(M_.params(11));
oo_.steady_state( 6 ) = log(M_.params(13));
oo_.steady_state( 3 ) = log(M_.params(10));
oo_.steady_state( 7 ) = log(M_.params(14));
oo_.steady_state( 5 ) = log(M_.params(12));
oo_.steady_state( 1 ) = log(M_.params(8));
oo_.steady_state( 2 ) = log(M_.params(9));
oo_.steady_state( 8 ) = log(M_.params(15));
oo_.steady_state( 9 ) = M_.params(16);
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.order = 1;
var_list_=[];
var_list_ = 'Y';
var_list_ = char(var_list_, 'C');
var_list_ = char(var_list_, 'pi');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'r');
var_list_ = char(var_list_, 'm_growth_ann');
info = stoch_simul(var_list_);
save('indeterminacry_error_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('indeterminacry_error_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('indeterminacry_error_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('indeterminacry_error_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('indeterminacry_error_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
