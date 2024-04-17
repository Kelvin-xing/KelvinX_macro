%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'RBC_stylized_facts';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'RBC_stylized_facts.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'ea';
M_.exo_names_tex = 'ea';
M_.endo_names = 'a';
M_.endo_names_tex = 'a';
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names = char(M_.endo_names, 'yn');
M_.endo_names_tex = char(M_.endo_names_tex, 'yn');
M_.param_names = 'theta';
M_.param_names_tex = 'theta';
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'sda');
M_.param_names_tex = char(M_.param_names_tex, 'sda');
M_.param_names = char(M_.param_names, 'ns');
M_.param_names_tex = char(M_.param_names_tex, 'ns');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names = char(M_.param_names, 'ks');
M_.param_names_tex = char(M_.param_names_tex, 'ks');
M_.param_names = char(M_.param_names, 'as');
M_.param_names_tex = char(M_.param_names_tex, 'as');
M_.param_names = char(M_.param_names, 'is');
M_.param_names_tex = char(M_.param_names_tex, 'is');
M_.param_names = char(M_.param_names, 'rs');
M_.param_names_tex = char(M_.param_names_tex, 'rs');
M_.param_names = char(M_.param_names, 'ws');
M_.param_names_tex = char(M_.param_names_tex, 'ws');
M_.param_names = char(M_.param_names, 'Rs');
M_.param_names_tex = char(M_.param_names_tex, 'Rs');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 15;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('RBC_stylized_facts_dynamic');
M_.lead_lag_incidence = [
 1 3 0;
 0 4 0;
 0 5 13;
 2 6 0;
 0 7 0;
 0 8 0;
 0 9 0;
 0 10 14;
 0 11 0;
 0 12 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = repmat(NaN,15, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 33;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 2 ) = 0.33;
alpha = M_.params( 2 );
M_.params( 4 ) = 0.99;
beta = M_.params( 4 );
M_.params( 3 ) = 0.025;
delta = M_.params( 3 );
M_.params( 5 ) = 1;
sigma = M_.params( 5 );
M_.params( 6 ) = 0.974;
rho = M_.params( 6 );
M_.params( 7 ) = 0.009;
sda = M_.params( 7 );
M_.params( 11 ) = 1;
as = M_.params( 11 );
M_.params( 8 ) = 0.3333333333333333;
ns = M_.params( 8 );
M_.params( 13 ) = 1/M_.params(4);
rs = M_.params( 13 );
kn = (alpha/(1/beta - 1 +delta))^(1/(1-alpha));
M_.params( 10 ) = M_.params(8)*kn;
ks = M_.params( 10 );
M_.params( 15 ) = 1/M_.params(4)-1+M_.params(3);
Rs = M_.params( 15 );
M_.params( 14 ) = (1-M_.params(2))*kn^M_.params(2);
ws = M_.params( 14 );
ys =as*kn^alpha*ns;
M_.params( 12 ) = M_.params(3)*M_.params(10);
is = M_.params( 12 );
M_.params( 9 ) = ys-M_.params(12);
cs = M_.params( 9 );
M_.params( 1 ) = kn^M_.params(2)*(1-M_.params(2))*(1-M_.params(8))/M_.params(9)*M_.params(11);
theta = M_.params( 1 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = log(M_.params(10));
oo_.steady_state( 9 ) = log(ys);
oo_.steady_state( 3 ) = log(M_.params(9));
oo_.steady_state( 5 ) = log(M_.params(12));
oo_.steady_state( 1 ) = log(M_.params(11));
oo_.steady_state( 6 ) = log(M_.params(13));
oo_.steady_state( 8 ) = M_.params(15);
oo_.steady_state( 7 ) = log(M_.params(14));
oo_.steady_state( 2 ) = log(M_.params(8));
oo_.steady_state( 10 ) = log(ys/M_.params(8));
oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];
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
M_.Sigma_e(1, 1) = M_.params(7)^2;
M_.sigma_e_is_diagonal = 1;
resid(1);
steady;
check;
options_.hp_filter = 1600;
options_.order = 1;
options_.periods = 1000;
var_list_=[];
info = stoch_simul(var_list_);
save('RBC_stylized_facts_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
