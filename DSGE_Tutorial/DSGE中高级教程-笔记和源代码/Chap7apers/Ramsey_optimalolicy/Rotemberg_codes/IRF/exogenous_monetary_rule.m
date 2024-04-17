%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'exogenous_monetary_rule';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('exogenous_monetary_rule.log');
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
M_.endo_names = char(M_.endo_names, 'Util');
M_.endo_names_tex = char(M_.endo_names_tex, 'Util');
M_.endo_names_long = char(M_.endo_names_long, 'Util');
M_.endo_names = char(M_.endo_names, 'Welf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Welf');
M_.endo_names_long = char(M_.endo_names_long, 'Welf');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.param_names = 'nbeta';
M_.param_names_tex = 'nbeta';
M_.param_names_long = 'nbeta';
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
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
M_.param_names = char(M_.param_names, 'A_SS');
M_.param_names_tex = char(M_.param_names_tex, 'A\_SS');
M_.param_names_long = char(M_.param_names_long, 'A_SS');
M_.param_names = char(M_.param_names, 'h_SS');
M_.param_names_tex = char(M_.param_names_tex, 'h\_SS');
M_.param_names_long = char(M_.param_names_long, 'h_SS');
M_.param_names = char(M_.param_names, 'R_SS');
M_.param_names_tex = char(M_.param_names_tex, 'R\_SS');
M_.param_names_long = char(M_.param_names_long, 'R_SS');
M_.param_names = char(M_.param_names, 'pie_SS');
M_.param_names_tex = char(M_.param_names_tex, 'pie\_SS');
M_.param_names_long = char(M_.param_names_long, 'pie_SS');
M_.param_names = char(M_.param_names, 'C_SS');
M_.param_names_tex = char(M_.param_names_tex, 'C\_SS');
M_.param_names_long = char(M_.param_names_long, 'C_SS');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 7;
M_.param_nbr = 14;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('exogenous_monetary_rule_static');
erase_compiled_function('exogenous_monetary_rule_dynamic');
M_.lead_lag_incidence = [
 0 2 0;
 0 3 0;
 0 4 9;
 0 5 10;
 0 6 0;
 0 7 11;
 1 8 0;]';
M_.nstatic = 3;
M_.nfwrd   = 3;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 1;
M_.ndynamic   = 4;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(14, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 24;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 3 ) = 0.99;
beta = M_.params( 3 );
M_.params( 4 ) = 5;
epsilon = M_.params( 4 );
M_.params( 5 ) = 100;
phi = M_.params( 5 );
M_.params( 6 ) = 0.9;
rho = M_.params( 6 );
M_.params( 7 ) = 1.5;
alpha = M_.params( 7 );
M_.params( 1 ) = M_.params(3);
nbeta = M_.params( 1 );
load parameterfile_irf;
set_param_value('nu',nu);
M_.params( 2 ) = 1;
chi = M_.params( 2 );
M_.params( 9 ) = 1;
pietarget = M_.params( 9 );
M_.params( 10 ) = 1;
A_SS = M_.params( 10 );
M_.params( 12 ) = M_.params(9)/M_.params(3)-1;
R_SS = M_.params( 12 );
M_.params( 13 ) = M_.params(9);
pie_SS = M_.params( 13 );
M_.params( 11 ) = ((1+M_.params(8))*(M_.params(4)-1)+M_.params(9)*M_.params(5)*(M_.params(9)-1)*(1-M_.params(3))/(1+M_.params(5)/2*(M_.params(9)-1)^2))/(M_.params(4)*M_.params(2));
h_SS = M_.params( 11 );
M_.params( 11 ) = sqrt(M_.params(11));
h_SS = M_.params( 11 );
M_.params( 14 ) = M_.params(11)/(1+M_.params(5)/2*(M_.params(13)-1)^2);
C_SS = M_.params( 14 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = M_.params(12);
oo_.steady_state( 2 ) = M_.params(11);
oo_.steady_state( 3 ) = M_.params(13);
oo_.steady_state( 4 ) = M_.params(14);
oo_.steady_state( 7 ) = M_.params(10);
oo_.steady_state( 5 ) = log(oo_.steady_state(4))-M_.params(2)*oo_.steady_state(2)^2/2;
oo_.steady_state( 6 ) = oo_.steady_state(5)/(1-M_.params(1));
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
options_.irf = 40;
options_.nograph = 1;
options_.order = 1;
var_list_=[];
var_list_ = 'pie';
var_list_ = char(var_list_, 'h');
var_list_ = char(var_list_, 'R');
var_list_ = char(var_list_, 'C');
var_list_ = char(var_list_, 'A');
var_list_ = char(var_list_, 'Util');
var_list_ = char(var_list_, 'Welf');
info = stoch_simul(var_list_);
save('exogenous_monetary_rule_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('exogenous_monetary_rule_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('exogenous_monetary_rule_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('exogenous_monetary_rule_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('exogenous_monetary_rule_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
