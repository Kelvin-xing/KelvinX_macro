%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'BGG1';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'BGG1.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'e_a';
M_.exo_names_tex = 'e\_a';
M_.exo_names = char(M_.exo_names, 'e_g');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_g');
M_.exo_names = char(M_.exo_names, 'e_rn');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_rn');
M_.endo_names = 'cH';
M_.endo_names_tex = 'cH';
M_.endo_names = char(M_.endo_names, 'hH');
M_.endo_names_tex = char(M_.endo_names_tex, 'hH');
M_.endo_names = char(M_.endo_names, 'piH');
M_.endo_names_tex = char(M_.endo_names_tex, 'piH');
M_.endo_names = char(M_.endo_names, 'rH');
M_.endo_names_tex = char(M_.endo_names_tex, 'rH');
M_.endo_names = char(M_.endo_names, 'r_nH');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_nH');
M_.endo_names = char(M_.endo_names, 'qH');
M_.endo_names_tex = char(M_.endo_names_tex, 'qH');
M_.endo_names = char(M_.endo_names, 'kH');
M_.endo_names_tex = char(M_.endo_names_tex, 'kH');
M_.endo_names = char(M_.endo_names, 'nH');
M_.endo_names_tex = char(M_.endo_names_tex, 'nH');
M_.endo_names = char(M_.endo_names, 'r_kH');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_kH');
M_.endo_names = char(M_.endo_names, 'yH');
M_.endo_names_tex = char(M_.endo_names_tex, 'yH');
M_.endo_names = char(M_.endo_names, 'xH');
M_.endo_names_tex = char(M_.endo_names_tex, 'xH');
M_.endo_names = char(M_.endo_names, 'iH');
M_.endo_names_tex = char(M_.endo_names_tex, 'iH');
M_.endo_names = char(M_.endo_names, 'aH');
M_.endo_names_tex = char(M_.endo_names_tex, 'aH');
M_.endo_names = char(M_.endo_names, 'c_eH');
M_.endo_names_tex = char(M_.endo_names_tex, 'c\_eH');
M_.endo_names = char(M_.endo_names, 'gH');
M_.endo_names_tex = char(M_.endo_names_tex, 'gH');
M_.endo_names = char(M_.endo_names, 'pi_t1H');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_t1H');
M_.endo_names = char(M_.endo_names, 'premiumH');
M_.endo_names_tex = char(M_.endo_names_tex, 'premiumH');
M_.endo_names = char(M_.endo_names, 'AUX_ENDO_LEAD_82');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_ENDO\_LEAD\_82');
M_.param_names = 'X';
M_.param_names_tex = 'X';
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names = char(M_.param_names, 'H');
M_.param_names_tex = char(M_.param_names_tex, 'H');
M_.param_names = char(M_.param_names, 'R_K');
M_.param_names_tex = char(M_.param_names_tex, 'R\_K');
M_.param_names = char(M_.param_names, 's');
M_.param_names_tex = char(M_.param_names_tex, 's');
M_.param_names = char(M_.param_names, 'KN');
M_.param_names_tex = char(M_.param_names_tex, 'KN');
M_.param_names = char(M_.param_names, 'CY');
M_.param_names_tex = char(M_.param_names_tex, 'CY');
M_.param_names = char(M_.param_names, 'GY');
M_.param_names_tex = char(M_.param_names_tex, 'GY');
M_.param_names = char(M_.param_names, 'C_EY');
M_.param_names_tex = char(M_.param_names_tex, 'C\_EY');
M_.param_names = char(M_.param_names, 'IY');
M_.param_names_tex = char(M_.param_names_tex, 'IY');
M_.param_names = char(M_.param_names, 'YK');
M_.param_names_tex = char(M_.param_names_tex, 'YK');
M_.param_names = char(M_.param_names, 'WY');
M_.param_names_tex = char(M_.param_names_tex, 'WY');
M_.param_names = char(M_.param_names, 'GAMMA_WBAR');
M_.param_names_tex = char(M_.param_names_tex, 'GAMMA\_WBAR');
M_.param_names = char(M_.param_names, 'NY');
M_.param_names_tex = char(M_.param_names_tex, 'NY');
M_.param_names = char(M_.param_names, 'DY');
M_.param_names_tex = char(M_.param_names_tex, 'DY');
M_.param_names = char(M_.param_names, 'YN');
M_.param_names_tex = char(M_.param_names_tex, 'YN');
M_.param_names = char(M_.param_names, 'niv');
M_.param_names_tex = char(M_.param_names_tex, 'niv');
M_.param_names = char(M_.param_names, 'omegav');
M_.param_names_tex = char(M_.param_names_tex, 'omegav');
M_.param_names = char(M_.param_names, 'alphav');
M_.param_names_tex = char(M_.param_names_tex, 'alphav');
M_.param_names = char(M_.param_names, 'betav');
M_.param_names_tex = char(M_.param_names_tex, 'betav');
M_.param_names = char(M_.param_names, 'sigmav');
M_.param_names_tex = char(M_.param_names_tex, 'sigmav');
M_.param_names = char(M_.param_names, 'gammav');
M_.param_names_tex = char(M_.param_names_tex, 'gammav');
M_.param_names = char(M_.param_names, 'muv');
M_.param_names_tex = char(M_.param_names_tex, 'muv');
M_.param_names = char(M_.param_names, 'deltav');
M_.param_names_tex = char(M_.param_names_tex, 'deltav');
M_.param_names = char(M_.param_names, 'phiv');
M_.param_names_tex = char(M_.param_names_tex, 'phiv');
M_.param_names = char(M_.param_names, 'kappav');
M_.param_names_tex = char(M_.param_names_tex, 'kappav');
M_.param_names = char(M_.param_names, 'thetav');
M_.param_names_tex = char(M_.param_names_tex, 'thetav');
M_.param_names = char(M_.param_names, 'epsilonv');
M_.param_names_tex = char(M_.param_names_tex, 'epsilonv');
M_.param_names = char(M_.param_names, 'etav');
M_.param_names_tex = char(M_.param_names_tex, 'etav');
M_.param_names = char(M_.param_names, 'zetav');
M_.param_names_tex = char(M_.param_names_tex, 'zetav');
M_.param_names = char(M_.param_names, 'rhov');
M_.param_names_tex = char(M_.param_names_tex, 'rhov');
M_.param_names = char(M_.param_names, 'rhov_a');
M_.param_names_tex = char(M_.param_names_tex, 'rhov\_a');
M_.param_names = char(M_.param_names, 'rhov_g');
M_.param_names_tex = char(M_.param_names_tex, 'rhov\_g');
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 18;
M_.param_nbr = 33;
M_.orig_endo_nbr = 17;
M_.aux_vars(1).endo_index = 18;
M_.aux_vars(1).type = 0;
M_.Sigma_e = zeros(3, 3);
M_.H = 0;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('BGG1_static');
erase_compiled_function('BGG1_dynamic');
M_.lead_lag_incidence = [
 0 9 27;
 0 10 0;
 0 11 28;
 1 12 0;
 2 13 0;
 3 14 0;
 4 15 0;
 5 16 0;
 0 17 29;
 0 18 0;
 0 19 30;
 0 20 0;
 6 21 0;
 0 22 0;
 7 23 0;
 8 24 0;
 0 25 0;
 0 26 31;]';
M_.nstatic = 5;
M_.nfwrd   = 5;
M_.npred   = 8;
M_.nboth   = 0;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(18, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(33, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 66;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
M_.params( 19 ) = 0.35;
alphav = M_.params( 19 );
M_.params( 18 ) = 0.64/(1-M_.params(19));
omegav = M_.params( 18 );
M_.params( 20 ) = 0.99;
betav = M_.params( 20 );
M_.params( 21 ) = 0.28;
sigmav = M_.params( 21 );
M_.params( 22 ) = 0.9728;
gammav = M_.params( 22 );
M_.params( 23 ) = 0.12;
muv = M_.params( 23 );
M_.params( 24 ) = 0.025;
deltav = M_.params( 24 );
M_.params( 31 ) = 0.9;
rhov = M_.params( 31 );
M_.params( 32 ) = 1;
rhov_a = M_.params( 32 );
M_.params( 33 ) = 0.95;
rhov_g = M_.params( 33 );
M_.params( 5 ) = 1.005;
s = M_.params( 5 );
M_.params( 6 ) = 2;
KN = M_.params( 6 );
M_.params( 8 ) = 0.2;
GY = M_.params( 8 );
M_.params( 1 ) = 1.1;
X = M_.params( 1 );
M_.params( 3 ) = 0.25;
H = M_.params( 3 );
M_.params( 2 ) = 1/M_.params(20);
R = M_.params( 2 );
M_.params( 4 ) = M_.params(5)*M_.params(2);
R_K = M_.params( 4 );
M_.params( 11 ) = M_.params(1)/M_.params(19)*(M_.params(4)-(1-M_.params(24)));
YK = M_.params( 11 );
M_.params( 12 ) = (1-M_.params(19))*M_.params(18)/(M_.params(1)*M_.params(3));
WY = M_.params( 12 );
M_.params( 13 ) = 1-1/(M_.params(5)*M_.params(22))*(M_.params(20)/M_.params(6)-(1-M_.params(19))*(1-M_.params(18))/M_.params(19)*(M_.params(5)-(1-M_.params(20)*M_.params(24))));
GAMMA_WBAR = M_.params( 13 );
M_.params( 14 ) = M_.params(4)*M_.params(22)*(1-M_.params(13))/M_.params(11)+(1-M_.params(19))*(1-M_.params(18))/M_.params(1);
NY = M_.params( 14 );
M_.params( 15 ) = 1/M_.params(11)-M_.params(14);
DY = M_.params( 15 );
M_.params( 7 ) = M_.params(3)*M_.params(12)-M_.params(8)+(M_.params(1)-1)/M_.params(1)+(M_.params(2)-1)*M_.params(15);
CY = M_.params( 7 );
M_.params( 10 ) = M_.params(24)*1/M_.params(11);
IY = M_.params( 10 );
M_.params( 9 ) = 1-M_.params(7)-M_.params(10)-M_.params(8);
C_EY = M_.params( 9 );
M_.params( 16 ) = 1/M_.params(14);
YN = M_.params( 16 );
M_.params( 17 ) = 0.05;
niv = M_.params( 17 );
M_.params( 28 ) = (1-M_.params(24))/(1-M_.params(24)+M_.params(19)*M_.params(11)/M_.params(1));
epsilonv = M_.params( 28 );
M_.params( 25 ) = 0.25;
phiv = M_.params( 25 );
M_.params( 29 ) = 3;
etav = M_.params( 29 );
M_.params( 27 ) = 0.75;
thetav = M_.params( 27 );
M_.params( 26 ) = (1-M_.params(27))/M_.params(27)*(1-M_.params(20)*M_.params(27));
kappav = M_.params( 26 );
M_.params( 30 ) = 0.11;
zetav = M_.params( 30 );
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.0065)^2;
M_.Sigma_e(2, 2) = (0.01)^2;
M_.Sigma_e(3, 3) = (0.000625)^2;
M_.sigma_e_is_diagonal = 1;
options_.irf = 30;
options_.nograph = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('BGG1_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
diary off
