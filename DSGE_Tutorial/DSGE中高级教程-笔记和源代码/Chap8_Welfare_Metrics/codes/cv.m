%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'cv';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'cv.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.endo_names = 'w';
M_.endo_names_tex = 'w';
M_.endo_names = char(M_.endo_names, 'w_c');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_c');
M_.endo_names = char(M_.endo_names, 'w_l');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_l');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names = char(M_.endo_names, 'wage');
M_.endo_names_tex = char(M_.endo_names_tex, 'wage');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names = char(M_.endo_names, 'Rk');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rk');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'sigmae');
M_.param_names_tex = char(M_.param_names_tex, 'sigmae');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names = char(M_.param_names, 'ws');
M_.param_names_tex = char(M_.param_names_tex, 'ws');
M_.param_names = char(M_.param_names, 'wcs');
M_.param_names_tex = char(M_.param_names_tex, 'wcs');
M_.param_names = char(M_.param_names, 'wls');
M_.param_names_tex = char(M_.param_names_tex, 'wls');
M_.param_names = char(M_.param_names, 'ks');
M_.param_names_tex = char(M_.param_names_tex, 'ks');
M_.param_names = char(M_.param_names, 'as');
M_.param_names_tex = char(M_.param_names_tex, 'as');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names = char(M_.param_names, 'ns');
M_.param_names_tex = char(M_.param_names_tex, 'ns');
M_.param_names = char(M_.param_names, 'wages');
M_.param_names_tex = char(M_.param_names_tex, 'wages');
M_.param_names = char(M_.param_names, 'ys');
M_.param_names_tex = char(M_.param_names_tex, 'ys');
M_.param_names = char(M_.param_names, 'is');
M_.param_names_tex = char(M_.param_names_tex, 'is');
M_.param_names = char(M_.param_names, 'Rks');
M_.param_names_tex = char(M_.param_names_tex, 'Rks');
M_.param_names = char(M_.param_names, 'rs');
M_.param_names_tex = char(M_.param_names_tex, 'rs');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 12;
M_.param_nbr = 20;
M_.orig_endo_nbr = 12;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('cv_dynamic');
M_.lead_lag_incidence = [
 0 3 0;
 0 4 15;
 0 5 16;
 1 6 0;
 2 7 0;
 0 8 17;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 18;
 0 14 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(12, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = repmat(NaN,20, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 38;
M_.NNZDerivatives(2) = 34;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.3333333333333333;
alpha = M_.params( 1 );
M_.params( 2 ) = 0.995;
beta = M_.params( 2 );
M_.params( 3 ) = 0.02;
delta = M_.params( 3 );
M_.params( 4 ) = 0.95;
rho = M_.params( 4 );
M_.params( 6 ) = 1.05;
sigma = M_.params( 6 );
M_.params( 7 ) = 0.4;
phi = M_.params( 7 );
M_.params( 15 ) = 0.3333333333333333;
ns = M_.params( 15 );
M_.params( 12 ) = M_.params(15)*(M_.params(1)/(1/M_.params(2)-1+M_.params(3)))^(1/(1-M_.params(1)));
ks = M_.params( 12 );
M_.params( 19 ) = M_.params(1)*M_.params(12)^(M_.params(1)-1)*M_.params(15)^(1-M_.params(1));
Rks = M_.params( 19 );
M_.params( 20 ) = M_.params(19)-M_.params(3);
rs = M_.params( 20 );
M_.params( 16 ) = (1-M_.params(1))*M_.params(15)^(-M_.params(1))*M_.params(12)^M_.params(1);
wages = M_.params( 16 );
M_.params( 17 ) = M_.params(15)^(1-M_.params(1))*M_.params(12)^M_.params(1);
ys = M_.params( 17 );
M_.params( 14 ) = M_.params(17)-M_.params(3)*M_.params(12);
cs = M_.params( 14 );
M_.params( 18 ) = M_.params(3)*M_.params(12);
is = M_.params( 18 );
M_.params( 13 ) = 1;
as = M_.params( 13 );
M_.params( 8 ) = M_.params(16)*M_.params(14)^(-M_.params(6))/M_.params(15)^M_.params(7);
psi = M_.params( 8 );
M_.params( 10 ) = 1/(1-M_.params(2))*(M_.params(14)^(1-M_.params(6))-1)/(1-M_.params(6));
wcs = M_.params( 10 );
M_.params( 11 ) = (-1)/(1-M_.params(2))*M_.params(8)*M_.params(15)^(1+M_.params(7))/(1+M_.params(7));
wls = M_.params( 11 );
M_.params( 9 ) = M_.params(10)+M_.params(11);
ws = M_.params( 9 );
load parameterfile_cv;
set_param_value('sigmae',sigmae);
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = M_.params(12);
oo_.steady_state( 7 ) = M_.params(15);
oo_.steady_state( 11 ) = M_.params(19);
oo_.steady_state( 12 ) = M_.params(20);
oo_.steady_state( 8 ) = M_.params(16);
oo_.steady_state( 9 ) = M_.params(17);
oo_.steady_state( 6 ) = M_.params(14);
oo_.steady_state( 10 ) = M_.params(18);
oo_.steady_state( 5 ) = M_.params(13);
oo_.steady_state( 2 ) = M_.params(10);
oo_.steady_state( 3 ) = M_.params(11);
oo_.steady_state( 1 ) = M_.params(10)+M_.params(11);
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
M_.Sigma_e(1, 1) = M_.params(5)^2;
M_.sigma_e_is_diagonal = 1;
steady;
resid(1);
options_.irf = 0;
options_.noprint = 1;
options_.order = 2;
var_list_=[];
info = stoch_simul(var_list_);
save('cv_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
