%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'Rotemberg_DIY';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Rotemberg_DIY.log');
M_.exo_names = 'eps_A';
M_.exo_names_tex = 'eps\_A';
M_.exo_names_long = 'eps_A';
M_.endo_names = 'R';
M_.endo_names_tex = 'R';
M_.endo_names_long = 'R';
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'pie');
M_.endo_names_tex = char(M_.endo_names_tex, 'pie');
M_.endo_names_long = char(M_.endo_names_long, 'pie');
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, 'C');
M_.endo_names_long = char(M_.endo_names_long, 'C');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.endo_names = char(M_.endo_names, 'Util');
M_.endo_names_tex = char(M_.endo_names_tex, 'Util');
M_.endo_names_long = char(M_.endo_names_long, 'Util');
M_.endo_names = char(M_.endo_names, 'Welf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Welf');
M_.endo_names_long = char(M_.endo_names_long, 'Welf');
M_.endo_names = char(M_.endo_names, 'lmult1');
M_.endo_names_tex = char(M_.endo_names_tex, 'lmult1');
M_.endo_names_long = char(M_.endo_names_long, 'lmult1');
M_.endo_names = char(M_.endo_names, 'lmult2');
M_.endo_names_tex = char(M_.endo_names_tex, 'lmult2');
M_.endo_names_long = char(M_.endo_names_long, 'lmult2');
M_.endo_names = char(M_.endo_names, 'lmult3');
M_.endo_names_tex = char(M_.endo_names_tex, 'lmult3');
M_.endo_names_long = char(M_.endo_names_long, 'lmult3');
M_.param_names = 'nbeta';
M_.param_names_tex = 'nbeta';
M_.param_names_long = 'nbeta';
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'epsil');
M_.param_names_tex = char(M_.param_names_tex, 'epsil');
M_.param_names_long = char(M_.param_names_long, 'epsil');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'nu');
M_.param_names_tex = char(M_.param_names_tex, 'nu');
M_.param_names_long = char(M_.param_names_long, 'nu');
M_.param_names = char(M_.param_names, 'pietarget');
M_.param_names_tex = char(M_.param_names_tex, 'pietarget');
M_.param_names_long = char(M_.param_names_long, 'pietarget');
M_.param_names = char(M_.param_names, 'as');
M_.param_names_tex = char(M_.param_names_tex, 'as');
M_.param_names_long = char(M_.param_names_long, 'as');
M_.param_names = char(M_.param_names, 'rs');
M_.param_names_tex = char(M_.param_names_tex, 'rs');
M_.param_names_long = char(M_.param_names_long, 'rs');
M_.param_names = char(M_.param_names, 'pis');
M_.param_names_tex = char(M_.param_names_tex, 'pis');
M_.param_names_long = char(M_.param_names_long, 'pis');
M_.param_names = char(M_.param_names, 'hs');
M_.param_names_tex = char(M_.param_names_tex, 'hs');
M_.param_names_long = char(M_.param_names_long, 'hs');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names_long = char(M_.param_names_long, 'cs');
M_.param_names = char(M_.param_names, 'Utils');
M_.param_names_tex = char(M_.param_names_tex, 'Utils');
M_.param_names_long = char(M_.param_names_long, 'Utils');
M_.param_names = char(M_.param_names, 'lmult1_SS');
M_.param_names_tex = char(M_.param_names_tex, 'lmult1\_SS');
M_.param_names_long = char(M_.param_names_long, 'lmult1_SS');
M_.param_names = char(M_.param_names, 'lmult2_SS');
M_.param_names_tex = char(M_.param_names_tex, 'lmult2\_SS');
M_.param_names_long = char(M_.param_names_long, 'lmult2_SS');
M_.param_names = char(M_.param_names, 'lmult3_SS');
M_.param_names_tex = char(M_.param_names_tex, 'lmult3\_SS');
M_.param_names_long = char(M_.param_names_long, 'lmult3_SS');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 18;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('Rotemberg_DIY_static');
erase_compiled_function('Rotemberg_DIY_dynamic');
M_.lead_lag_incidence = [
 0 5 0;
 0 6 0;
 0 7 15;
 1 8 16;
 2 9 0;
 0 10 0;
 0 11 17;
 3 12 0;
 4 13 0;
 0 14 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(18, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 50;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 3 ) = 0.99;
beta = M_.params( 3 );
M_.params( 4 ) = 5;
epsil = M_.params( 4 );
M_.params( 5 ) = 100;
phi = M_.params( 5 );
M_.params( 6 ) = 0.9;
rho = M_.params( 6 );
M_.params( 7 ) = 1.5;
alpha = M_.params( 7 );
M_.params( 1 ) = M_.params(3);
nbeta = M_.params( 1 );
M_.params( 8 ) = 0;
nu = M_.params( 8 );
M_.params( 2 ) = 1;
chi = M_.params( 2 );
M_.params( 9 ) = 1.;
pietarget = M_.params( 9 );
M_.params( 10 ) = 1;
as = M_.params( 10 );
M_.params( 11 ) = M_.params(9)/M_.params(3)-1;
rs = M_.params( 11 );
M_.params( 12 ) = M_.params(9);
pis = M_.params( 12 );
M_.params( 13 ) = ((1+M_.params(8))*(M_.params(4)-1)+M_.params(9)*M_.params(5)*(M_.params(9)-1)*(1-M_.params(3))/(1+M_.params(5)/2*(M_.params(9)-1)^2))/(M_.params(4)*M_.params(2));
hs = M_.params( 13 );
M_.params( 13 ) = sqrt(M_.params(13));
hs = M_.params( 13 );
M_.params( 14 ) = M_.params(13)/(1+M_.params(5)/2*(M_.params(12)-1)^2);
cs = M_.params( 14 );
M_.params( 15 ) = log(M_.params(14))-M_.params(2)*M_.params(13)^2/2;
Utils = M_.params( 15 );
M_.params( 16 ) = 0;
lmult1_SS = M_.params( 16 );
M_.params( 17 ) = 0;
lmult2_SS = M_.params( 17 );
M_.params( 18 ) = (-1);
lmult3_SS = M_.params( 18 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = M_.params(11);
oo_.steady_state( 2 ) = M_.params(13);
oo_.steady_state( 3 ) = M_.params(12);
oo_.steady_state( 4 ) = M_.params(14);
oo_.steady_state( 5 ) = 1;
oo_.steady_state( 8 ) = M_.params(16);
oo_.steady_state( 9 ) = M_.params(17);
oo_.steady_state( 10 ) = M_.params(18);
oo_.steady_state( 6 ) = log(M_.params(14))-M_.params(2)*M_.params(13)^2/2;
oo_.steady_state( 7 ) = log(M_.params(14))-M_.params(2)*M_.params(13)^2/2/(1-M_.params(1));
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
M_.Sigma_e(1, 1) = (.01)^2;
options_.irf = 20;
options_.order = 1;
options_.periods = 200;
var_list_=[];
var_list_ = 'pie';
var_list_ = char(var_list_, 'h');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'C');
var_list_ = char(var_list_, 'A');
var_list_ = char(var_list_, 'lmult1');
var_list_ = char(var_list_, 'lmult2');
var_list_ = char(var_list_, 'lmult3');
var_list_ = char(var_list_, 'Welf');
var_list_ = char(var_list_, 'Util');
info = stoch_simul(var_list_);
save('Rotemberg_DIY_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Rotemberg_DIY_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Rotemberg_DIY_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Rotemberg_DIY_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Rotemberg_DIY_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
