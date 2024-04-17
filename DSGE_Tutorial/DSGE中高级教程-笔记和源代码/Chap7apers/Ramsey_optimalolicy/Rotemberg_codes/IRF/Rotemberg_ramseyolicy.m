%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'Rotemberg_ramsey_policy';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Rotemberg_ramsey_policy.log');
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
M_.endo_names = char(M_.endo_names, 'MULT_1');
M_.endo_names_tex = char(M_.endo_names_tex, 'MULT\_1');
M_.endo_names_long = char(M_.endo_names_long, 'MULT_1');
M_.endo_names = char(M_.endo_names, 'MULT_2');
M_.endo_names_tex = char(M_.endo_names_tex, 'MULT\_2');
M_.endo_names_long = char(M_.endo_names_long, 'MULT_2');
M_.endo_names = char(M_.endo_names, 'MULT_3');
M_.endo_names_tex = char(M_.endo_names_tex, 'MULT\_3');
M_.endo_names_long = char(M_.endo_names_long, 'MULT_3');
M_.endo_names = char(M_.endo_names, 'MULT_4');
M_.endo_names_tex = char(M_.endo_names_tex, 'MULT\_4');
M_.endo_names_long = char(M_.endo_names_long, 'MULT_4');
M_.endo_names = char(M_.endo_names, 'MULT_5');
M_.endo_names_tex = char(M_.endo_names_tex, 'MULT\_5');
M_.endo_names_long = char(M_.endo_names_long, 'MULT_5');
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
M_.param_names = char(M_.param_names, 'optimal_policy_discount_factor');
M_.param_names_tex = char(M_.param_names_tex, 'optimal\_policy\_discount\_factor');
M_.param_names_long = char(M_.param_names_long, 'optimal_policy_discount_factor');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 11;
M_.param_nbr = 15;
M_.orig_endo_nbr = 6;
M_.aux_vars(1).endo_index = 7;
M_.aux_vars(1).type = 6;
M_.aux_vars(1).eq_nbr = '1';
M_.aux_vars(2).endo_index = 8;
M_.aux_vars(2).type = 6;
M_.aux_vars(2).eq_nbr = '2';
M_.aux_vars(3).endo_index = 9;
M_.aux_vars(3).type = 6;
M_.aux_vars(3).eq_nbr = '3';
M_.aux_vars(4).endo_index = 10;
M_.aux_vars(4).type = 6;
M_.aux_vars(4).eq_nbr = '4';
M_.aux_vars(5).endo_index = 11;
M_.aux_vars(5).type = 6;
M_.aux_vars(5).eq_nbr = '5';
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('Rotemberg_ramsey_policy_static');
erase_compiled_function('Rotemberg_ramsey_policy_dynamic');
M_.orig_eq_nbr = 5;
M_.lead_lag_incidence = [
 0 5 0;
 0 6 0;
 0 7 16;
 1 8 17;
 2 9 0;
 0 10 0;
 3 11 0;
 4 12 0;
 0 13 0;
 0 14 18;
 0 15 0;]';
M_.nstatic = 5;
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
oo_.steady_state = zeros(11, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(15, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 58;
M_.NNZDerivatives(2) = 160;
M_.NNZDerivatives(3) = -1;
M_.params( 3 ) = 0.99;
beta = M_.params( 3 );
M_.params( 4 ) = 5;
epsilon = M_.params( 4 );
M_.params( 5 ) = 100;
phi = M_.params( 5 );
M_.params( 6 ) = 0.9;
rho = M_.params( 6 );
M_.params( 1 ) = M_.params(3);
nbeta = M_.params( 1 );
load parameterfile_irf;
set_param_value('nu',nu);
M_.params( 2 ) = 1;
chi = M_.params( 2 );
M_.params( 8 ) = 1.;
pietarget = M_.params( 8 );
M_.params( 9 ) = 1;
as = M_.params( 9 );
M_.params( 10 ) = M_.params(8)/M_.params(3)-1;
rs = M_.params( 10 );
M_.params( 11 ) = M_.params(8);
pis = M_.params( 11 );
M_.params( 12 ) = ((1+M_.params(7))*(M_.params(4)-1)+M_.params(11)*M_.params(5)*(M_.params(11)-1)*(1-M_.params(3))/(1+M_.params(5)/2*(M_.params(11)-1)^2))/(M_.params(4)*M_.params(2));
hs = M_.params( 12 );
M_.params( 12 ) = sqrt(M_.params(12));
hs = M_.params( 12 );
M_.params( 13 ) = M_.params(12)/(1+M_.params(5)/2*(M_.params(11)-1)^2);
cs = M_.params( 13 );
M_.params( 14 ) = log(M_.params(13))-M_.params(2)*M_.params(12)^2/2;
Utils = M_.params( 14 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 6 ) = M_.params(14);
oo_.steady_state( 1 ) = M_.params(10);
oo_.steady_state( 2 ) = M_.params(12);
oo_.steady_state( 3 ) = M_.params(11);
oo_.steady_state( 4 ) = M_.params(13);
oo_.steady_state( 5 ) = M_.params(9);
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
M_.params( 15 ) = 0.99;
optimal_policy_discount_factor = M_.params( 15 );
options_.irf = 40;
options_.nograph = 1;
options_.order = 1;
var_list_=[];
ramsey_policy(var_list_);
save('Rotemberg_ramsey_policy_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Rotemberg_ramsey_policy_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Rotemberg_ramsey_policy_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Rotemberg_ramsey_policy_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Rotemberg_ramsey_policy_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
