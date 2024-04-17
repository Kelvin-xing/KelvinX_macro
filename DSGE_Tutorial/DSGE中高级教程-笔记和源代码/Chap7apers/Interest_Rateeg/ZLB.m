%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'ZLB';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'ea';
M_.exo_names_tex = 'ea';
M_.exo_names_long = 'ea';
M_.exo_names = char(M_.exo_names, 'ei');
M_.exo_names_tex = char(M_.exo_names_tex, 'ei');
M_.exo_names_long = char(M_.exo_names_long, 'ei');
M_.exo_names = char(M_.exo_names, 'eg');
M_.exo_names_tex = char(M_.exo_names_tex, 'eg');
M_.exo_names_long = char(M_.exo_names_long, 'eg');
M_.exo_names = char(M_.exo_names, 'ex');
M_.exo_names_tex = char(M_.exo_names_tex, 'ex');
M_.exo_names_long = char(M_.exo_names_long, 'ex');
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'yf');
M_.endo_names_tex = char(M_.endo_names_tex, 'yf');
M_.endo_names_long = char(M_.endo_names_long, 'yf');
M_.endo_names = char(M_.endo_names, 'vp');
M_.endo_names_tex = char(M_.endo_names_tex, 'vp');
M_.endo_names_long = char(M_.endo_names_long, 'vp');
M_.endo_names = char(M_.endo_names, 'pisharp');
M_.endo_names_tex = char(M_.endo_names_tex, 'pisharp');
M_.endo_names_long = char(M_.endo_names_long, 'pisharp');
M_.endo_names = char(M_.endo_names, 'x1');
M_.endo_names_tex = char(M_.endo_names_tex, 'x1');
M_.endo_names_long = char(M_.endo_names_long, 'x1');
M_.endo_names = char(M_.endo_names, 'x2');
M_.endo_names_tex = char(M_.endo_names_tex, 'x2');
M_.endo_names_long = char(M_.endo_names_long, 'x2');
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'g');
M_.endo_names = char(M_.endo_names, 'nex');
M_.endo_names_tex = char(M_.endo_names_tex, 'nex');
M_.endo_names_long = char(M_.endo_names_long, 'nex');
M_.param_names = 'sigma';
M_.param_names_tex = 'sigma';
M_.param_names_long = 'sigma';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'phipi');
M_.param_names_tex = char(M_.param_names_tex, 'phipi');
M_.param_names_long = char(M_.param_names_long, 'phipi');
M_.param_names = char(M_.param_names, 'phiy');
M_.param_names_tex = char(M_.param_names_tex, 'phiy');
M_.param_names_long = char(M_.param_names_long, 'phiy');
M_.param_names = char(M_.param_names, 'etag');
M_.param_names_tex = char(M_.param_names_tex, 'etag');
M_.param_names_long = char(M_.param_names_long, 'etag');
M_.param_names = char(M_.param_names, 'etax');
M_.param_names_tex = char(M_.param_names_tex, 'etax');
M_.param_names_long = char(M_.param_names_long, 'etax');
M_.param_names = char(M_.param_names, 'rhog');
M_.param_names_tex = char(M_.param_names_tex, 'rhog');
M_.param_names_long = char(M_.param_names_long, 'rhog');
M_.param_names = char(M_.param_names, 'rhox');
M_.param_names_tex = char(M_.param_names_tex, 'rhox');
M_.param_names_long = char(M_.param_names_long, 'rhox');
M_.param_names = char(M_.param_names, 'rhoa');
M_.param_names_tex = char(M_.param_names_tex, 'rhoa');
M_.param_names_long = char(M_.param_names_long, 'rhoa');
M_.param_names = char(M_.param_names, 'rhoi');
M_.param_names_tex = char(M_.param_names_tex, 'rhoi');
M_.param_names_long = char(M_.param_names_long, 'rhoi');
M_.param_names = char(M_.param_names, 'sigmaa');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaa');
M_.param_names_long = char(M_.param_names_long, 'sigmaa');
M_.param_names = char(M_.param_names, 'sigmai');
M_.param_names_tex = char(M_.param_names_tex, 'sigmai');
M_.param_names_long = char(M_.param_names_long, 'sigmai');
M_.param_names = char(M_.param_names, 'sigmax');
M_.param_names_tex = char(M_.param_names_tex, 'sigmax');
M_.param_names_long = char(M_.param_names_long, 'sigmax');
M_.param_names = char(M_.param_names, 'sigmag');
M_.param_names_tex = char(M_.param_names_tex, 'sigmag');
M_.param_names_long = char(M_.param_names_long, 'sigmag');
M_.param_names = char(M_.param_names, 'cs');
M_.param_names_tex = char(M_.param_names_tex, 'cs');
M_.param_names_long = char(M_.param_names_long, 'cs');
M_.param_names = char(M_.param_names, 'is');
M_.param_names_tex = char(M_.param_names_tex, 'is');
M_.param_names_long = char(M_.param_names_long, 'is');
M_.param_names = char(M_.param_names, 'pis');
M_.param_names_tex = char(M_.param_names_tex, 'pis');
M_.param_names_long = char(M_.param_names_long, 'pis');
M_.param_names = char(M_.param_names, 'ns');
M_.param_names_tex = char(M_.param_names_tex, 'ns');
M_.param_names_long = char(M_.param_names_long, 'ns');
M_.param_names = char(M_.param_names, 'ws');
M_.param_names_tex = char(M_.param_names_tex, 'ws');
M_.param_names_long = char(M_.param_names_long, 'ws');
M_.param_names = char(M_.param_names, 'mcs');
M_.param_names_tex = char(M_.param_names_tex, 'mcs');
M_.param_names_long = char(M_.param_names_long, 'mcs');
M_.param_names = char(M_.param_names, 'as');
M_.param_names_tex = char(M_.param_names_tex, 'as');
M_.param_names_long = char(M_.param_names_long, 'as');
M_.param_names = char(M_.param_names, 'ys');
M_.param_names_tex = char(M_.param_names_tex, 'ys');
M_.param_names_long = char(M_.param_names_long, 'ys');
M_.param_names = char(M_.param_names, 'yfs');
M_.param_names_tex = char(M_.param_names_tex, 'yfs');
M_.param_names_long = char(M_.param_names_long, 'yfs');
M_.param_names = char(M_.param_names, 'vps');
M_.param_names_tex = char(M_.param_names_tex, 'vps');
M_.param_names_long = char(M_.param_names_long, 'vps');
M_.param_names = char(M_.param_names, 'pisharps');
M_.param_names_tex = char(M_.param_names_tex, 'pisharps');
M_.param_names_long = char(M_.param_names_long, 'pisharps');
M_.param_names = char(M_.param_names, 'x1s');
M_.param_names_tex = char(M_.param_names_tex, 'x1s');
M_.param_names_long = char(M_.param_names_long, 'x1s');
M_.param_names = char(M_.param_names, 'x2s');
M_.param_names_tex = char(M_.param_names_tex, 'x2s');
M_.param_names_long = char(M_.param_names_long, 'x2s');
M_.param_names = char(M_.param_names, 'gs');
M_.param_names_tex = char(M_.param_names_tex, 'gs');
M_.param_names_long = char(M_.param_names_long, 'gs');
M_.param_names = char(M_.param_names, 'nexs');
M_.param_names_tex = char(M_.param_names_tex, 'nexs');
M_.param_names_long = char(M_.param_names_long, 'nexs');
M_.exo_det_nbr = 0;
M_.exo_nbr = 4;
M_.endo_nbr = 16;
M_.param_nbr = 33;
M_.orig_endo_nbr = 16;
M_.aux_vars = [];
M_.Sigma_e = zeros(4, 4);
M_.Correlation_matrix = eye(4, 4);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('ZLB_static');
erase_compiled_function('ZLB_dynamic');
M_.lead_lag_incidence = [
 0 6 22;
 1 7 0;
 0 8 23;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 2 13 0;
 0 14 0;
 0 15 0;
 3 16 0;
 0 17 0;
 0 18 24;
 0 19 25;
 4 20 0;
 5 21 0;]';
M_.nstatic = 7;
M_.nfwrd   = 4;
M_.npred   = 5;
M_.nboth   = 0;
M_.nsfwrd   = 4;
M_.nspred   = 5;
M_.ndynamic   = 9;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:4];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(16, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(4, 1);
M_.params = NaN(33, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 59;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 2 ) = .99;
beta = M_.params( 2 );
M_.params( 1 ) = 1;
sigma = M_.params( 1 );
M_.params( 4 ) = 1;
eta = M_.params( 4 );
M_.params( 3 ) = 1;
psi = M_.params( 3 );
M_.params( 6 ) = 10;
epsilon = M_.params( 6 );
M_.params( 7 ) = 1.5;
phipi = M_.params( 7 );
M_.params( 8 ) = 0;
phiy = M_.params( 8 );
M_.params( 9 ) = .1413;
etag = M_.params( 9 );
M_.params( 10 ) = .0320;
etax = M_.params( 10 );
M_.params( 13 ) = .95;
rhoa = M_.params( 13 );
M_.params( 14 ) = .0;
rhoi = M_.params( 14 );
M_.params( 11 ) = .9658;
rhog = M_.params( 11 );
M_.params( 12 ) = .8257;
rhox = M_.params( 12 );
M_.params( 15 ) = .1;
sigmaa = M_.params( 15 );
M_.params( 16 ) = .1;
sigmai = M_.params( 16 );
M_.params( 17 ) = .5010;
sigmax = M_.params( 17 );
M_.params( 18 ) = .1459;
sigmag = M_.params( 18 );
M_.params( 5 ) = .75;
phi = M_.params( 5 );
M_.params( 21 ) = 1;
pis = M_.params( 21 );
M_.params( 20 ) = 1/M_.params(2)*M_.params(21);
is = M_.params( 20 );
M_.params( 25 ) = 1;
as = M_.params( 25 );
M_.params( 29 ) = ((M_.params(21)^(1-M_.params(6))-M_.params(5))/(1-M_.params(5)))^(1/(1-M_.params(6)));
pisharps = M_.params( 29 );
M_.params( 28 ) = (1-M_.params(5))*(M_.params(21)/M_.params(29))^M_.params(6)/(1-M_.params(5)*M_.params(21)^M_.params(6));
vps = M_.params( 28 );
M_.params( 24 ) = (M_.params(6)-1)*M_.params(29)*(1-M_.params(21)^M_.params(6)*M_.params(2)*M_.params(5))/(1-M_.params(2)*M_.params(5)*M_.params(21)^(M_.params(6)-1))/M_.params(21)/M_.params(6);
mcs = M_.params( 24 );
M_.params( 22 ) = (M_.params(28)^M_.params(1)*M_.params(24)/M_.params(3)/(1-M_.params(10)-M_.params(9))^M_.params(1))^(1/(M_.params(1)+M_.params(4)));
ns = M_.params( 22 );
M_.params( 26 ) = M_.params(25)*M_.params(22)/M_.params(28);
ys = M_.params( 26 );
M_.params( 27 ) = ((M_.params(6)-1)/M_.params(6)/M_.params(3)/(1-M_.params(9)-M_.params(10))^M_.params(1))^(1/(M_.params(1)+M_.params(4)))*M_.params(25)^((1+M_.params(4))/(M_.params(1)+M_.params(4)));
yfs = M_.params( 27 );
M_.params( 32 ) = M_.params(9)*M_.params(26);
gs = M_.params( 32 );
M_.params( 33 ) = M_.params(10)*M_.params(26);
nexs = M_.params( 33 );
M_.params( 19 ) = (1-M_.params(9)-M_.params(10))*M_.params(26);
cs = M_.params( 19 );
M_.params( 23 ) = M_.params(24)*M_.params(25);
ws = M_.params( 23 );
M_.params( 30 ) = M_.params(26)*M_.params(24)*M_.params(19)^(-M_.params(1))/(1-M_.params(21)^M_.params(6)*M_.params(2)*M_.params(5));
x1s = M_.params( 30 );
M_.params( 31 ) = M_.params(26)*M_.params(19)^(-M_.params(1))/(1-M_.params(2)*M_.params(5)*M_.params(21)^(M_.params(6)-1));
x2s = M_.params( 31 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = log(M_.params(19));
oo_.steady_state( 2 ) = log(M_.params(20));
oo_.steady_state( 3 ) = log(M_.params(21));
oo_.steady_state( 5 ) = log(M_.params(22));
oo_.steady_state( 6 ) = log(M_.params(23));
oo_.steady_state( 7 ) = log(M_.params(24));
oo_.steady_state( 8 ) = log(M_.params(25));
oo_.steady_state( 9 ) = log(M_.params(26));
oo_.steady_state( 10 ) = log(M_.params(27));
oo_.steady_state( 11 ) = log(M_.params(28));
oo_.steady_state( 12 ) = log(M_.params(29));
oo_.steady_state( 13 ) = log(M_.params(30));
oo_.steady_state( 14 ) = log(M_.params(31));
oo_.steady_state( 15 ) = log(M_.params(32));
oo_.steady_state( 4 ) = log(M_.params(20)/M_.params(21));
oo_.steady_state( 16 ) = log(M_.params(33));
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
M_.Sigma_e(1, 1) = M_.params(15)^2;
M_.Sigma_e(2, 2) = M_.params(16)^2;
M_.Sigma_e(3, 3) = M_.params(18)^2;
M_.Sigma_e(4, 4) = M_.params(17)^2;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.nograph = 1;
options_.noprint = 1;
options_.order = 1;
var_list_=[];
var_list_ = 'i';
var_list_ = char(var_list_, 'pi');
var_list_ = char(var_list_, 'n');
var_list_ = char(var_list_, 'w');
var_list_ = char(var_list_, 'mc');
var_list_ = char(var_list_, 'y');
var_list_ = char(var_list_, 'r');
var_list_ = char(var_list_, 'g');
var_list_ = char(var_list_, 'c');
var_list_ = char(var_list_, 'nex');
var_list_ = char(var_list_, 'yf');
info = stoch_simul(var_list_);
save('ZLB_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('ZLB_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('ZLB_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('ZLB_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('ZLB_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
