%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'MIU';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'MIU.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'ea';
M_.exo_names_tex = 'ea';
M_.exo_names = char(M_.exo_names, 'em');
M_.exo_names_tex = char(M_.exo_names_tex, 'em');
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names = char(M_.endo_names, 'dlnm');
M_.endo_names_tex = char(M_.endo_names_tex, 'dlnm');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.param_names = 'theta';
M_.param_names_tex = 'theta';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names = char(M_.param_names, 'rhom');
M_.param_names_tex = char(M_.param_names_tex, 'rhom');
M_.param_names = char(M_.param_names, 'pistar');
M_.param_names_tex = char(M_.param_names_tex, 'pistar');
M_.param_names = char(M_.param_names, 'sigmam');
M_.param_names_tex = char(M_.param_names_tex, 'sigmam');
M_.param_names = char(M_.param_names, 'sigmaa');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaa');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names = char(M_.param_names, 'rhoa');
M_.param_names_tex = char(M_.param_names_tex, 'rhoa');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names = char(M_.param_names, 'ns');
M_.param_names_tex = char(M_.param_names_tex, 'ns');
M_.param_names = char(M_.param_names, 'ws');
M_.param_names_tex = char(M_.param_names_tex, 'ws');
M_.param_names = char(M_.param_names, 'Rs');
M_.param_names_tex = char(M_.param_names_tex, 'Rs');
M_.param_names = char(M_.param_names, 'is');
M_.param_names_tex = char(M_.param_names_tex, 'is');
M_.param_names = char(M_.param_names, 'ps');
M_.param_names_tex = char(M_.param_names_tex, 'ps');
M_.param_names = char(M_.param_names, 'ks');
M_.param_names_tex = char(M_.param_names_tex, 'ks');
M_.param_names = char(M_.param_names, 'ys');
M_.param_names_tex = char(M_.param_names_tex, 'ys');
M_.param_names = char(M_.param_names, 'Is');
M_.param_names_tex = char(M_.param_names_tex, 'Is');
M_.param_names = char(M_.param_names, 'ms');
M_.param_names_tex = char(M_.param_names_tex, 'ms');
M_.param_names = char(M_.param_names, 'rs');
M_.param_names_tex = char(M_.param_names_tex, 'rs');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 13;
M_.param_nbr = 22;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('MIU_dynamic');
M_.lead_lag_incidence = [
 0 6 19;
 0 7 0;
 0 8 0;
 0 9 20;
 0 10 0;
 1 11 0;
 2 12 0;
 0 13 0;
 0 14 0;
 3 15 0;
 4 16 21;
 5 17 0;
 0 18 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = repmat(NaN,22, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 45;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 2 ) = .99;
beta = M_.params( 2 );
M_.params( 3 ) = 0.3333333333333333;
alpha = M_.params( 3 );
M_.params( 4 ) = .025;
delta = M_.params( 4 );
M_.params( 5 ) = 1;
psi = M_.params( 5 );
M_.params( 10 ) = 1;
zeta = M_.params( 10 );
M_.params( 6 ) = .5;
rhom = M_.params( 6 );
M_.params( 11 ) = .5;
rhoa = M_.params( 11 );
M_.params( 8 ) = .01;
sigmam = M_.params( 8 );
M_.params( 9 ) = .01;
sigmaa = M_.params( 9 );
M_.params( 7 ) = 1.02;
pistar = M_.params( 7 );
M_.params( 15 ) = 1/M_.params(2)-1+M_.params(4);
Rs = M_.params( 15 );
kn = (alpha/Rs)^(1/(1-alpha));
M_.params( 14 ) = (1-M_.params(3))*kn^M_.params(3);
ws = M_.params( 14 );
M_.params( 16 ) = M_.params(7)/M_.params(2);
is = M_.params( 16 );
M_.params( 13 ) = 0.3333333333333333;
ns = M_.params( 13 );
M_.params( 18 ) = kn*M_.params(13);
ks = M_.params( 18 );
M_.params( 20 ) = M_.params(4)*M_.params(18);
Is = M_.params( 20 );
M_.params( 19 ) = kn^M_.params(3)*M_.params(13);
ys = M_.params( 19 );
M_.params( 12 ) = M_.params(19)-M_.params(20);
cs = M_.params( 12 );
M_.params( 21 ) = M_.params(5)^M_.params(10)*M_.params(12)^M_.params(10)*(M_.params(16)/(M_.params(16)-1))^M_.params(10);
ms = M_.params( 21 );
M_.params( 22 ) = M_.params(16)/M_.params(7);
rs = M_.params( 22 );
M_.params( 1 ) = M_.params(14)/M_.params(12)*(1-M_.params(13));
theta = M_.params( 1 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = log(M_.params(12));
oo_.steady_state( 2 ) = log(M_.params(13));
oo_.steady_state( 3 ) = log(M_.params(14));
oo_.steady_state( 4 ) = log(M_.params(15));
oo_.steady_state( 5 ) = log(M_.params(16));
oo_.steady_state( 6 ) = log(M_.params(18));
oo_.steady_state( 7 ) = 0;
oo_.steady_state( 8 ) = log(M_.params(19));
oo_.steady_state( 9 ) = log(M_.params(20));
oo_.steady_state( 10 ) = 0;
oo_.steady_state( 11 ) = log(M_.params(7));
oo_.steady_state( 12 ) = log(M_.params(21));
oo_.steady_state( 13 ) = log(M_.params(22));
oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
resid;
steady;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(9)^2;
M_.Sigma_e(2, 2) = M_.params(8)^2;
M_.sigma_e_is_diagonal = 1;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('MIU_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
