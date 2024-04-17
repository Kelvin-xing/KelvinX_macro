%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'recursively_running_dynare';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('recursively_running_dynare.log');
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.exo_names_long = 'e';
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names_long = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names_long = char(M_.param_names_long, 'sigma');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'sigmae');
M_.param_names_tex = char(M_.param_names_tex, 'sigmae');
M_.param_names_long = char(M_.param_names_long, 'sigmae');
M_.param_names = char(M_.param_names, 'yss');
M_.param_names_tex = char(M_.param_names_tex, 'yss');
M_.param_names_long = char(M_.param_names_long, 'yss');
M_.param_names = char(M_.param_names, 'iss');
M_.param_names_tex = char(M_.param_names_tex, 'iss');
M_.param_names_long = char(M_.param_names_long, 'iss');
M_.param_names = char(M_.param_names, 'rss');
M_.param_names_tex = char(M_.param_names_tex, 'rss');
M_.param_names_long = char(M_.param_names_long, 'rss');
M_.param_names = char(M_.param_names, 'nss');
M_.param_names_tex = char(M_.param_names_tex, 'nss');
M_.param_names_long = char(M_.param_names_long, 'nss');
M_.param_names = char(M_.param_names, 'wss');
M_.param_names_tex = char(M_.param_names_tex, 'wss');
M_.param_names_long = char(M_.param_names_long, 'wss');
M_.param_names = char(M_.param_names, 'kss');
M_.param_names_tex = char(M_.param_names_tex, 'kss');
M_.param_names_long = char(M_.param_names_long, 'kss');
M_.param_names = char(M_.param_names, 'Rss');
M_.param_names_tex = char(M_.param_names_tex, 'Rss');
M_.param_names_long = char(M_.param_names_long, 'Rss');
M_.param_names = char(M_.param_names, 'css');
M_.param_names_tex = char(M_.param_names_tex, 'css');
M_.param_names_long = char(M_.param_names_long, 'css');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 9;
M_.param_nbr = 16;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('recursively_running_dynare_static');
erase_compiled_function('recursively_running_dynare_dynamic');
M_.lead_lag_incidence = [
 0 3 0;
 0 4 0;
 0 5 0;
 0 6 0;
 0 7 0;
 1 8 0;
 2 9 0;
 0 10 12;
 0 11 13;]';
M_.nstatic = 5;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(16, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 30;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
load parametersaved;
set_param_value('alpha',alpha);
M_.params( 2 ) = .99;
beta = M_.params( 2 );
M_.params( 3 ) = .02;
delta = M_.params( 3 );
M_.params( 4 ) = .97;
rho = M_.params( 4 );
M_.params( 5 ) = 1;
sigma = M_.params( 5 );
M_.params( 6 ) = 3;
psi = M_.params( 6 );
M_.params( 7 ) = 1;
eta = M_.params( 7 );
M_.params( 8 ) = .01;
sigmae = M_.params( 8 );
M_.params( 11 ) = 1/M_.params(2)-1;
rss = M_.params( 11 );
M_.params( 15 ) = 1/M_.params(2)-1+M_.params(3);
Rss = M_.params( 15 );
kn = (alpha/Rss)^(1/(1-alpha));
M_.params( 13 ) = (1-M_.params(1))*kn^M_.params(1);
wss = M_.params( 13 );
M_.params( 12 ) = ((kn^M_.params(1)-M_.params(3)*kn)*(M_.params(6)/M_.params(13))^(1/M_.params(5)))^((-1)/(1+M_.params(7)/M_.params(5)));
nss = M_.params( 12 );
M_.params( 14 ) = kn*M_.params(12);
kss = M_.params( 14 );
M_.params( 9 ) = kn^M_.params(1)*M_.params(12);
yss = M_.params( 9 );
M_.params( 10 ) = M_.params(3)*M_.params(14);
iss = M_.params( 10 );
M_.params( 16 ) = M_.params(9)-M_.params(10);
css = M_.params( 16 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 6 ) = log(M_.params(14));
oo_.steady_state( 1 ) = log(M_.params(9));
oo_.steady_state( 9 ) = log(M_.params(16));
oo_.steady_state( 2 ) = log(M_.params(10));
oo_.steady_state( 7 ) = 0;
oo_.steady_state( 3 ) = M_.params(11);
oo_.steady_state( 8 ) = M_.params(15);
oo_.steady_state( 5 ) = log(M_.params(13));
oo_.steady_state( 4 ) = log(M_.params(12));
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
M_.Sigma_e(1, 1) = M_.params(8)^2;
resid(1); 
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.nograph = 1;
options_.noprint = 1;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('recursively_running_dynare_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('recursively_running_dynare_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('recursively_running_dynare_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('recursively_running_dynare_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('recursively_running_dynare_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
