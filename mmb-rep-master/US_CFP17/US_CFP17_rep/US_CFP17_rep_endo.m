%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'US_CFP17_rep_endo';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('US_CFP17_rep_endo.log');
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names_long = 'eps_a';
M_.exo_names = char(M_.exo_names, 'eps_mp');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_mp');
M_.exo_names_long = char(M_.exo_names_long, 'eps_mp');
M_.exo_names = char(M_.exo_names, 'eps_i');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_i');
M_.exo_names_long = char(M_.exo_names_long, 'eps_i');
M_.exo_names = char(M_.exo_names, 'eps_psi');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_psi');
M_.exo_names_long = char(M_.exo_names_long, 'eps_psi');
M_.exo_names = char(M_.exo_names, 'eps_mk');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_mk');
M_.exo_names_long = char(M_.exo_names_long, 'eps_mk');
M_.exo_names = char(M_.exo_names, 'eps_mkw');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_mkw');
M_.exo_names_long = char(M_.exo_names_long, 'eps_mkw');
M_.exo_names = char(M_.exo_names, 'eps_b2');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_b2');
M_.exo_names_long = char(M_.exo_names_long, 'eps_b2');
M_.exo_names = char(M_.exo_names, 'eps_rn');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_rn');
M_.exo_names_long = char(M_.exo_names_long, 'eps_rn');
M_.endo_names = 'rn';
M_.endo_names_tex = 'rn';
M_.endo_names_long = 'rn';
M_.endo_names = char(M_.endo_names, 'mrs');
M_.endo_names_tex = char(M_.endo_names_tex, 'mrs');
M_.endo_names_long = char(M_.endo_names_long, 'mrs');
M_.endo_names = char(M_.endo_names, 'pinw');
M_.endo_names_tex = char(M_.endo_names_tex, 'pinw');
M_.endo_names_long = char(M_.endo_names_long, 'pinw');
M_.endo_names = char(M_.endo_names, 'mkw');
M_.endo_names_tex = char(M_.endo_names_tex, 'mkw');
M_.endo_names_long = char(M_.endo_names_long, 'mkw');
M_.endo_names = char(M_.endo_names, 'mk');
M_.endo_names_tex = char(M_.endo_names_tex, 'mk');
M_.endo_names_long = char(M_.endo_names_long, 'mk');
M_.endo_names = char(M_.endo_names, 'gy');
M_.endo_names_tex = char(M_.endo_names_tex, 'gy');
M_.endo_names_long = char(M_.endo_names_long, 'gy');
M_.endo_names = char(M_.endo_names, 'gi');
M_.endo_names_tex = char(M_.endo_names_tex, 'gi');
M_.endo_names_long = char(M_.endo_names_long, 'gi');
M_.endo_names = char(M_.endo_names, 'spread');
M_.endo_names_tex = char(M_.endo_names_tex, 'spread');
M_.endo_names_long = char(M_.endo_names_long, 'spread');
M_.endo_names = char(M_.endo_names, 'g_bonds');
M_.endo_names_tex = char(M_.endo_names_tex, 'g\_bonds');
M_.endo_names_long = char(M_.endo_names_long, 'g_bonds');
M_.endo_names = char(M_.endo_names, 'infl');
M_.endo_names_tex = char(M_.endo_names_tex, 'infl');
M_.endo_names_long = char(M_.endo_names_long, 'infl');
M_.endo_names = char(M_.endo_names, 'ffr');
M_.endo_names_tex = char(M_.endo_names_tex, 'ffr');
M_.endo_names_long = char(M_.endo_names_long, 'ffr');
M_.endo_names = char(M_.endo_names, 'bb2');
M_.endo_names_tex = char(M_.endo_names_tex, 'bb2');
M_.endo_names_long = char(M_.endo_names_long, 'bb2');
M_.endo_names = char(M_.endo_names, 'qnat');
M_.endo_names_tex = char(M_.endo_names_tex, 'qnat');
M_.endo_names_long = char(M_.endo_names_long, 'qnat');
M_.endo_names = char(M_.endo_names, 'r10s');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10s');
M_.endo_names_long = char(M_.endo_names_long, 'r10s');
M_.endo_names = char(M_.endo_names, 'assets');
M_.endo_names_tex = char(M_.endo_names_tex, 'assets');
M_.endo_names_long = char(M_.endo_names_long, 'assets');
M_.endo_names = char(M_.endo_names, 'qi');
M_.endo_names_tex = char(M_.endo_names_tex, 'qi');
M_.endo_names_long = char(M_.endo_names_long, 'qi');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names_long = char(M_.endo_names_long, 'm');
M_.endo_names = char(M_.endo_names, 'lev');
M_.endo_names_tex = char(M_.endo_names_tex, 'lev');
M_.endo_names_long = char(M_.endo_names_long, 'lev');
M_.endo_names = char(M_.endo_names, 'muc');
M_.endo_names_tex = char(M_.endo_names_tex, 'muc');
M_.endo_names_long = char(M_.endo_names_long, 'muc');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'r1');
M_.endo_names_tex = char(M_.endo_names_tex, 'r1');
M_.endo_names_long = char(M_.endo_names_long, 'r1');
M_.endo_names = char(M_.endo_names, 'r2');
M_.endo_names_tex = char(M_.endo_names_tex, 'r2');
M_.endo_names_long = char(M_.endo_names_long, 'r2');
M_.endo_names = char(M_.endo_names, 'rk');
M_.endo_names_tex = char(M_.endo_names_tex, 'rk');
M_.endo_names_long = char(M_.endo_names_long, 'rk');
M_.endo_names = char(M_.endo_names, 'pin');
M_.endo_names_tex = char(M_.endo_names_tex, 'pin');
M_.endo_names_long = char(M_.endo_names_long, 'pin');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'mp');
M_.endo_names_tex = char(M_.endo_names_tex, 'mp');
M_.endo_names_long = char(M_.endo_names_long, 'mp');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'f');
M_.endo_names_tex = char(M_.endo_names_tex, 'f');
M_.endo_names_long = char(M_.endo_names_long, 'f');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'b2');
M_.endo_names_tex = char(M_.endo_names_tex, 'b2');
M_.endo_names_long = char(M_.endo_names_long, 'b2');
M_.endo_names = char(M_.endo_names, 'nw');
M_.endo_names_tex = char(M_.endo_names_tex, 'nw');
M_.endo_names_long = char(M_.endo_names_long, 'nw');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'pk');
M_.endo_names_tex = char(M_.endo_names_tex, 'pk');
M_.endo_names_long = char(M_.endo_names_long, 'pk');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'ann_pin');
M_.endo_names_tex = char(M_.endo_names_tex, 'ann\_pin');
M_.endo_names_long = char(M_.endo_names_long, 'ann_pin');
M_.endo_names = char(M_.endo_names, 'r10');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10');
M_.endo_names_long = char(M_.endo_names_long, 'r10');
M_.endo_names = char(M_.endo_names, 'ann_r2');
M_.endo_names_tex = char(M_.endo_names_tex, 'ann\_r2');
M_.endo_names_long = char(M_.endo_names_long, 'ann_r2');
M_.endo_names = char(M_.endo_names, 'ann_r1');
M_.endo_names_tex = char(M_.endo_names_tex, 'ann\_r1');
M_.endo_names_long = char(M_.endo_names_long, 'ann_r1');
M_.endo_names = char(M_.endo_names, 'term_prem');
M_.endo_names_tex = char(M_.endo_names_tex, 'term\_prem');
M_.endo_names_long = char(M_.endo_names_long, 'term_prem');
M_.endo_names = char(M_.endo_names, 'r10_nat');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10\_nat');
M_.endo_names_long = char(M_.endo_names_long, 'r10_nat');
M_.endo_names = char(M_.endo_names, 'muinv');
M_.endo_names_tex = char(M_.endo_names_tex, 'muinv');
M_.endo_names_long = char(M_.endo_names_long, 'muinv');
M_.endo_names = char(M_.endo_names, 'r_lend');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_lend');
M_.endo_names_long = char(M_.endo_names_long, 'r_lend');
M_.endo_names = char(M_.endo_names, 'u_psi');
M_.endo_names_tex = char(M_.endo_names_tex, 'u\_psi');
M_.endo_names_long = char(M_.endo_names_long, 'u_psi');
M_.endo_names = char(M_.endo_names, 'mucf');
M_.endo_names_tex = char(M_.endo_names_tex, 'mucf');
M_.endo_names_long = char(M_.endo_names_long, 'mucf');
M_.endo_names = char(M_.endo_names, 'cf');
M_.endo_names_tex = char(M_.endo_names_tex, 'cf');
M_.endo_names_long = char(M_.endo_names_long, 'cf');
M_.endo_names = char(M_.endo_names, 'r1f');
M_.endo_names_tex = char(M_.endo_names_tex, 'r1f');
M_.endo_names_long = char(M_.endo_names_long, 'r1f');
M_.endo_names = char(M_.endo_names, 'pinf');
M_.endo_names_tex = char(M_.endo_names_tex, 'pinf');
M_.endo_names_long = char(M_.endo_names_long, 'pinf');
M_.endo_names = char(M_.endo_names, 'pkf');
M_.endo_names_tex = char(M_.endo_names_tex, 'pkf');
M_.endo_names_long = char(M_.endo_names_long, 'pkf');
M_.endo_names = char(M_.endo_names, 'mf');
M_.endo_names_tex = char(M_.endo_names_tex, 'mf');
M_.endo_names_long = char(M_.endo_names_long, 'mf');
M_.endo_names = char(M_.endo_names, 'rkf');
M_.endo_names_tex = char(M_.endo_names_tex, 'rkf');
M_.endo_names_long = char(M_.endo_names_long, 'rkf');
M_.endo_names = char(M_.endo_names, 'qif');
M_.endo_names_tex = char(M_.endo_names_tex, 'qif');
M_.endo_names_long = char(M_.endo_names_long, 'qif');
M_.endo_names = char(M_.endo_names, 'r2f');
M_.endo_names_tex = char(M_.endo_names_tex, 'r2f');
M_.endo_names_long = char(M_.endo_names_long, 'r2f');
M_.endo_names = char(M_.endo_names, 'Lf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lf');
M_.endo_names_long = char(M_.endo_names_long, 'Lf');
M_.endo_names = char(M_.endo_names, 'wf');
M_.endo_names_tex = char(M_.endo_names_tex, 'wf');
M_.endo_names_long = char(M_.endo_names_long, 'wf');
M_.endo_names = char(M_.endo_names, 'yf');
M_.endo_names_tex = char(M_.endo_names_tex, 'yf');
M_.endo_names_long = char(M_.endo_names_long, 'yf');
M_.endo_names = char(M_.endo_names, 'kf');
M_.endo_names_tex = char(M_.endo_names_tex, 'kf');
M_.endo_names_long = char(M_.endo_names_long, 'kf');
M_.endo_names = char(M_.endo_names, 'if');
M_.endo_names_tex = char(M_.endo_names_tex, 'if');
M_.endo_names_long = char(M_.endo_names_long, 'if');
M_.endo_names = char(M_.endo_names, 'ff');
M_.endo_names_tex = char(M_.endo_names_tex, 'ff');
M_.endo_names_long = char(M_.endo_names_long, 'ff');
M_.endo_names = char(M_.endo_names, 'nwf');
M_.endo_names_tex = char(M_.endo_names_tex, 'nwf');
M_.endo_names_long = char(M_.endo_names_long, 'nwf');
M_.endo_names = char(M_.endo_names, 'term_premf');
M_.endo_names_tex = char(M_.endo_names_tex, 'term\_premf');
M_.endo_names_long = char(M_.endo_names_long, 'term_premf');
M_.endo_names = char(M_.endo_names, 'r10f');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10f');
M_.endo_names_long = char(M_.endo_names_long, 'r10f');
M_.endo_names = char(M_.endo_names, 'r10_natf');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10\_natf');
M_.endo_names_long = char(M_.endo_names_long, 'r10_natf');
M_.endo_names = char(M_.endo_names, 'qnatf');
M_.endo_names_tex = char(M_.endo_names_tex, 'qnatf');
M_.endo_names_long = char(M_.endo_names_long, 'qnatf');
M_.endo_names = char(M_.endo_names, 'qf');
M_.endo_names_tex = char(M_.endo_names_tex, 'qf');
M_.endo_names_long = char(M_.endo_names_long, 'qf');
M_.endo_names = char(M_.endo_names, 'b2f');
M_.endo_names_tex = char(M_.endo_names_tex, 'b2f');
M_.endo_names_long = char(M_.endo_names_long, 'b2f');
M_.endo_names = char(M_.endo_names, 'ygap');
M_.endo_names_tex = char(M_.endo_names_tex, 'ygap');
M_.endo_names_long = char(M_.endo_names_long, 'ygap');
M_.endo_names = char(M_.endo_names, 'levf');
M_.endo_names_tex = char(M_.endo_names_tex, 'levf');
M_.endo_names_long = char(M_.endo_names_long, 'levf');
M_.endo_names = char(M_.endo_names, 'r10obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'r10obs');
M_.endo_names_long = char(M_.endo_names_long, 'r10obs');
M_.endo_names = char(M_.endo_names, 'PCE_inf');
M_.endo_names_tex = char(M_.endo_names_tex, 'PCE\_inf');
M_.endo_names_long = char(M_.endo_names_long, 'PCE_inf');
M_.endo_names = char(M_.endo_names, 'y_growth');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_growth');
M_.endo_names_long = char(M_.endo_names_long, 'y_growth');
M_.endo_names = char(M_.endo_names, 'i_growth');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_growth');
M_.endo_names_long = char(M_.endo_names_long, 'i_growth');
M_.endo_names = char(M_.endo_names, 'labor_dist');
M_.endo_names_tex = char(M_.endo_names_tex, 'labor\_dist');
M_.endo_names_long = char(M_.endo_names_long, 'labor_dist');
M_.endo_names = char(M_.endo_names, 'AUX_ENDO_LAG_11_1');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_ENDO\_LAG\_11\_1');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_ENDO_LAG_11_1');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names_long = 'alpha';
M_.param_names = char(M_.param_names, 'b');
M_.param_names_tex = char(M_.param_names_tex, 'b');
M_.param_names_long = char(M_.param_names_long, 'b');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'h');
M_.param_names_tex = char(M_.param_names_tex, 'h');
M_.param_names_long = char(M_.param_names_long, 'h');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'kappa_i');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_i');
M_.param_names_long = char(M_.param_names_long, 'kappa_i');
M_.param_names = char(M_.param_names, 'psi_i');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_i');
M_.param_names_long = char(M_.param_names_long, 'psi_i');
M_.param_names = char(M_.param_names, 'psi_n');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_n');
M_.param_names_long = char(M_.param_names_long, 'psi_n');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names_long = char(M_.param_names_long, 'zeta');
M_.param_names = char(M_.param_names, 'tau_p');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_p');
M_.param_names_long = char(M_.param_names_long, 'tau_p');
M_.param_names = char(M_.param_names, 'tau_y');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_y');
M_.param_names_long = char(M_.param_names_long, 'tau_y');
M_.param_names = char(M_.param_names, 'tau_pi');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_pi');
M_.param_names_long = char(M_.param_names_long, 'tau_pi');
M_.param_names = char(M_.param_names, 'tauy_long');
M_.param_names_tex = char(M_.param_names_tex, 'tauy\_long');
M_.param_names_long = char(M_.param_names_long, 'tauy_long');
M_.param_names = char(M_.param_names, 'taupi_long');
M_.param_names_tex = char(M_.param_names_tex, 'taupi\_long');
M_.param_names_long = char(M_.param_names_long, 'taupi_long');
M_.param_names = char(M_.param_names, 'tau_prem');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_prem');
M_.param_names_long = char(M_.param_names_long, 'tau_prem');
M_.param_names = char(M_.param_names, 'rho_m');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_m');
M_.param_names_long = char(M_.param_names_long, 'rho_m');
M_.param_names = char(M_.param_names, 'eps_p');
M_.param_names_tex = char(M_.param_names_tex, 'eps\_p');
M_.param_names_long = char(M_.param_names_long, 'eps_p');
M_.param_names = char(M_.param_names, 'eps_w');
M_.param_names_tex = char(M_.param_names_tex, 'eps\_w');
M_.param_names_long = char(M_.param_names_long, 'eps_w');
M_.param_names = char(M_.param_names, 'theta_p');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_p');
M_.param_names_long = char(M_.param_names_long, 'theta_p');
M_.param_names = char(M_.param_names, 'theta_w');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_w');
M_.param_names_long = char(M_.param_names_long, 'theta_w');
M_.param_names = char(M_.param_names, 'i_p');
M_.param_names_tex = char(M_.param_names_tex, 'i\_p');
M_.param_names_long = char(M_.param_names_long, 'i_p');
M_.param_names = char(M_.param_names, 'i_w');
M_.param_names_tex = char(M_.param_names_tex, 'i\_w');
M_.param_names_long = char(M_.param_names_long, 'i_w');
M_.param_names = char(M_.param_names, 'kappaw');
M_.param_names_tex = char(M_.param_names_tex, 'kappaw');
M_.param_names_long = char(M_.param_names_long, 'kappaw');
M_.param_names = char(M_.param_names, 'kappapc');
M_.param_names_tex = char(M_.param_names_tex, 'kappapc');
M_.param_names_long = char(M_.param_names_long, 'kappapc');
M_.param_names = char(M_.param_names, 'sigmamk');
M_.param_names_tex = char(M_.param_names_tex, 'sigmamk');
M_.param_names_long = char(M_.param_names_long, 'sigmamk');
M_.param_names = char(M_.param_names, 'sigmamkw');
M_.param_names_tex = char(M_.param_names_tex, 'sigmamkw');
M_.param_names_long = char(M_.param_names_long, 'sigmamkw');
M_.param_names = char(M_.param_names, 'sigmab2');
M_.param_names_tex = char(M_.param_names_tex, 'sigmab2');
M_.param_names_long = char(M_.param_names_long, 'sigmab2');
M_.param_names = char(M_.param_names, 'sigmarn');
M_.param_names_tex = char(M_.param_names_tex, 'sigmarn');
M_.param_names_long = char(M_.param_names_long, 'sigmarn');
M_.param_names = char(M_.param_names, 'sigmaea');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaea');
M_.param_names_long = char(M_.param_names_long, 'sigmaea');
M_.param_names = char(M_.param_names, 'sigmaeb');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaeb');
M_.param_names_long = char(M_.param_names_long, 'sigmaeb');
M_.param_names = char(M_.param_names, 'sigmaeph');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaeph');
M_.param_names_long = char(M_.param_names_long, 'sigmaeph');
M_.param_names = char(M_.param_names, 'sigmaer');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaer');
M_.param_names_long = char(M_.param_names_long, 'sigmaer');
M_.param_names = char(M_.param_names, 'sigmaemu');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaemu');
M_.param_names_long = char(M_.param_names_long, 'sigmaemu');
M_.param_names = char(M_.param_names, 'dur');
M_.param_names_tex = char(M_.param_names_tex, 'dur');
M_.param_names_long = char(M_.param_names_long, 'dur');
M_.param_names = char(M_.param_names, 'rhomkw');
M_.param_names_tex = char(M_.param_names_tex, 'rhomkw');
M_.param_names_long = char(M_.param_names_long, 'rhomkw');
M_.param_names = char(M_.param_names, 'rhomk');
M_.param_names_tex = char(M_.param_names_tex, 'rhomk');
M_.param_names_long = char(M_.param_names_long, 'rhomk');
M_.param_names = char(M_.param_names, 'rhod3');
M_.param_names_tex = char(M_.param_names_tex, 'rhod3');
M_.param_names_long = char(M_.param_names_long, 'rhod3');
M_.param_names = char(M_.param_names, 'rhod4');
M_.param_names_tex = char(M_.param_names_tex, 'rhod4');
M_.param_names_long = char(M_.param_names_long, 'rhod4');
M_.param_names = char(M_.param_names, 'rhoi');
M_.param_names_tex = char(M_.param_names_tex, 'rhoi');
M_.param_names_long = char(M_.param_names_long, 'rhoi');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'rho1_b');
M_.param_names_tex = char(M_.param_names_tex, 'rho1\_b');
M_.param_names_long = char(M_.param_names_long, 'rho1_b');
M_.param_names = char(M_.param_names, 'rho2_b');
M_.param_names_tex = char(M_.param_names_tex, 'rho2\_b');
M_.param_names_long = char(M_.param_names_long, 'rho2_b');
M_.param_names = char(M_.param_names, 'rho_phi');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_phi');
M_.param_names_long = char(M_.param_names_long, 'rho_phi');
M_.param_names = char(M_.param_names, 'rho_mu');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_mu');
M_.param_names_long = char(M_.param_names_long, 'rho_mu');
M_.param_names = char(M_.param_names, 'rho_rn');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_rn');
M_.param_names_long = char(M_.param_names_long, 'rho_rn');
M_.param_names = char(M_.param_names, 'rhoi_long');
M_.param_names_tex = char(M_.param_names_tex, 'rhoi\_long');
M_.param_names_long = char(M_.param_names_long, 'rhoi_long');
M_.param_names = char(M_.param_names, 'Y_ss');
M_.param_names_tex = char(M_.param_names_tex, 'Y\_ss');
M_.param_names_long = char(M_.param_names_long, 'Y_ss');
M_.param_names = char(M_.param_names, 'I_ss');
M_.param_names_tex = char(M_.param_names_tex, 'I\_ss');
M_.param_names_long = char(M_.param_names_long, 'I_ss');
M_.param_names = char(M_.param_names, 'C_ss');
M_.param_names_tex = char(M_.param_names_tex, 'C\_ss');
M_.param_names_long = char(M_.param_names_long, 'C_ss');
M_.param_names = char(M_.param_names, 'R1ss');
M_.param_names_tex = char(M_.param_names_tex, 'R1ss');
M_.param_names_long = char(M_.param_names_long, 'R1ss');
M_.param_names = char(M_.param_names, 'R2ss');
M_.param_names_tex = char(M_.param_names_tex, 'R2ss');
M_.param_names_long = char(M_.param_names_long, 'R2ss');
M_.param_names = char(M_.param_names, 'b2ss');
M_.param_names_tex = char(M_.param_names_tex, 'b2ss');
M_.param_names_long = char(M_.param_names_long, 'b2ss');
M_.param_names = char(M_.param_names, 'dss');
M_.param_names_tex = char(M_.param_names_tex, 'dss');
M_.param_names_long = char(M_.param_names_long, 'dss');
M_.param_names = char(M_.param_names, 'nwss');
M_.param_names_tex = char(M_.param_names_tex, 'nwss');
M_.param_names_long = char(M_.param_names_long, 'nwss');
M_.param_names = char(M_.param_names, 'xss');
M_.param_names_tex = char(M_.param_names_tex, 'xss');
M_.param_names_long = char(M_.param_names_long, 'xss');
M_.param_names = char(M_.param_names, 'premss');
M_.param_names_tex = char(M_.param_names_tex, 'premss');
M_.param_names_long = char(M_.param_names_long, 'premss');
M_.param_names = char(M_.param_names, 'b2n');
M_.param_names_tex = char(M_.param_names_tex, 'b2n');
M_.param_names_long = char(M_.param_names_long, 'b2n');
M_.exo_det_nbr = 0;
M_.exo_nbr = 8;
M_.endo_nbr = 77;
M_.param_nbr = 60;
M_.orig_endo_nbr = 76;
M_.aux_vars(1).endo_index = 77;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 12;
M_.aux_vars(1).orig_lead_lag = -1;
M_.Sigma_e = zeros(8, 8);
M_.Correlation_matrix = eye(8, 8);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('US_CFP17_rep_endo_static');
erase_compiled_function('US_CFP17_rep_endo_dynamic');
M_.lead_lag_incidence = [
 1 25 0;
 0 26 0;
 0 27 102;
 2 28 0;
 3 29 0;
 0 30 0;
 0 31 0;
 0 32 0;
 0 33 0;
 0 34 0;
 0 35 0;
 4 36 0;
 0 37 103;
 0 38 0;
 0 39 0;
 5 40 104;
 0 41 105;
 0 42 0;
 0 43 106;
 6 44 107;
 0 45 0;
 7 46 0;
 8 47 0;
 9 48 0;
 0 49 0;
 0 50 108;
 10 51 109;
 11 52 0;
 12 53 0;
 0 54 110;
 13 55 0;
 0 56 0;
 0 57 0;
 0 58 0;
 14 59 0;
 15 60 111;
 0 61 112;
 0 62 0;
 0 63 0;
 0 64 0;
 0 65 0;
 0 66 0;
 0 67 0;
 0 68 0;
 16 69 0;
 0 70 0;
 17 71 0;
 0 72 113;
 18 73 114;
 19 74 0;
 0 75 115;
 0 76 116;
 0 77 117;
 0 78 118;
 20 79 119;
 0 80 0;
 0 81 0;
 0 82 0;
 0 83 0;
 21 84 0;
 22 85 120;
 23 86 0;
 0 87 0;
 0 88 0;
 0 89 0;
 0 90 0;
 0 91 121;
 0 92 122;
 0 93 0;
 0 94 0;
 0 95 0;
 0 96 0;
 0 97 0;
 0 98 0;
 0 99 0;
 0 100 0;
 24 101 0;]';
M_.nstatic = 39;
M_.nfwrd   = 14;
M_.npred   = 17;
M_.nboth   = 7;
M_.nsfwrd   = 21;
M_.nspred   = 24;
M_.ndynamic   = 38;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:8];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(77, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(8, 1);
M_.params = NaN(60, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 272;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
load parameterfile;
set_param_value('alpha',alpha);           set_param_value('b',b);                   set_param_value('beta',beta);
set_param_value('delta',delta);           set_param_value('eta',eta);               set_param_value('gamma',gamma);
set_param_value('h',h);                   set_param_value('kappa',kappa);           set_param_value('kappa_i',kappa_i);
set_param_value('psi_i',psi_i);           set_param_value('psi_n',psi_n);           set_param_value('zeta',zeta);
set_param_value('tau_p',tau_p);           set_param_value('tau_y',tau_y);           set_param_value('tau_pi',tau_pi);
set_param_value('tau_prem',tau_prem);     set_param_value('tauy_long',tauy_long);   set_param_value('taupi_long',taupi_long);
set_param_value('nwss',nwss); 
set_param_value('theta_p',theta_p);       set_param_value('theta_w',theta_w);       set_param_value('i_p',i_p);
set_param_value('i_w',i_w);               set_param_value('kappaw',kappaw);         set_param_value('kappapc',kappapc);
set_param_value('sigmaea',sigmaea);       set_param_value('sigmaeb',sigmaeb);       set_param_value('sigmaeph',sigmaeph);
set_param_value('sigmaer',sigmaer);       set_param_value('sigmaemu',sigmaemu);     set_param_value('sigmarn',sigmarn);
set_param_value('sigmamk',sigmamk);       set_param_value('sigmamkw',sigmamkw);     set_param_value('sigmab2',sigmab2);
set_param_value('rhoi',rhoi);             set_param_value('rho_a',rho_a);           set_param_value('rho1_b',rho1_b);
set_param_value('rho2_b',rho2_b);         set_param_value('rho_phi',rho_phi);       set_param_value('rho_mu',rho_mu); set_param_value('rho_m',rho_m);
set_param_value('rho_rn',rho_rn);         set_param_value('rhomk',rhomk);           
set_param_value('rhomkw',rhomkw);         set_param_value('rhoi_long',rhoi_long);
set_param_value('C_ss',C_ss);             set_param_value('I_ss',I_ss);             set_param_value('R1ss',R1ss);                     
set_param_value('R2ss',R2ss);             set_param_value('Y_ss',Y_ss);             set_param_value('premss',premss); 
M_.params( 60 ) = .4*M_.params(6);
b2n = M_.params( 60 );
M_.params( 55 ) = M_.params(60)*M_.params(57);
b2ss = M_.params( 55 );
fn  = gamma - b2n;    fss  = nwss*fn;
M_.params( 58 ) = fss+M_.params(55);
xss = M_.params( 58 );
M_.params( 56 ) = M_.params(57)*(M_.params(6)-1);
dss = M_.params( 56 );
M_.params( 9 ) = (duration_i-1)/duration_i;
kappa_i = M_.params( 9 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 19 ) = 0;
oo_.steady_state( 8 ) = 0;
oo_.steady_state( 10 ) = 0;
oo_.steady_state( 11 ) = 0;
oo_.steady_state( 16 ) = 0;
oo_.steady_state( 18 ) = 0;
oo_.steady_state( 74 ) = 0;
oo_.steady_state( 73 ) = 0;
oo_.steady_state( 72 ) = 0;
oo_.steady_state( 27 ) = 0;
oo_.steady_state( 51 ) = 0;
oo_.steady_state( 36 ) = 0;
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 75 ) = 0;
oo_.steady_state( 5 ) = 0;
oo_.steady_state( 20 ) = 0;
oo_.steady_state( 21 ) = 0;
oo_.steady_state( 22 ) = 0;
oo_.steady_state( 23 ) = 0;
oo_.steady_state( 47 ) = 0;
oo_.steady_state( 45 ) = 0;
oo_.steady_state( 29 ) = 0;
oo_.steady_state( 40 ) = 0;
oo_.steady_state( 24 ) = 0;
oo_.steady_state( 25 ) = 0;
oo_.steady_state( 27 ) = 0;
oo_.steady_state( 28 ) = 0;
oo_.steady_state( 29 ) = 0;
oo_.steady_state( 30 ) = 0;
oo_.steady_state( 31 ) = 0;
oo_.steady_state( 32 ) = 0;
oo_.steady_state( 33 ) = 0;
oo_.steady_state( 34 ) = 0;
oo_.steady_state( 35 ) = 0;
oo_.steady_state( 36 ) = 0;
oo_.steady_state( 37 ) = 0;
oo_.steady_state( 38 ) = 0;
oo_.steady_state( 10 ) = 0;
oo_.steady_state( 11 ) = 0;
oo_.steady_state( 39 ) = 0;
oo_.steady_state( 40 ) = 0;
oo_.steady_state(77)=oo_.steady_state(12);
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(32)^2;
M_.Sigma_e(2, 2) = M_.params(35)^2;
M_.Sigma_e(3, 3) = M_.params(36)^2;
M_.Sigma_e(4, 4) = M_.params(34)^2;
M_.Sigma_e(5, 5) = M_.params(28)^2;
M_.Sigma_e(6, 6) = M_.params(29)^2;
M_.Sigma_e(7, 7) = M_.params(30)^2;
M_.Sigma_e(8, 8) = M_.params(31)^2;
options_.irf = 100;
options_.nograph = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('US_CFP17_rep_endo_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('US_CFP17_rep_endo_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('US_CFP17_rep_endo_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('US_CFP17_rep_endo_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('US_CFP17_rep_endo_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off