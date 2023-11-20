%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'US_PV15_rep';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('US_PV15_rep.log');
M_.parameter_used_with_lead_lag = true;
M_.exo_names = 'e_a';
M_.exo_names_tex = 'e\_a';
M_.exo_names_long = 'e_a';
M_.exo_names = char(M_.exo_names, 'e_g');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_g');
M_.exo_names_long = char(M_.exo_names_long, 'e_g');
M_.exo_names = char(M_.exo_names, 'e_b');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_b');
M_.exo_names_long = char(M_.exo_names_long, 'e_b');
M_.exo_names = char(M_.exo_names, 'e_i');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_i');
M_.exo_names_long = char(M_.exo_names_long, 'e_i');
M_.exo_names = char(M_.exo_names, 'e_l');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_l');
M_.exo_names_long = char(M_.exo_names_long, 'e_l');
M_.exo_names = char(M_.exo_names, 'e_n');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_n');
M_.exo_names_long = char(M_.exo_names_long, 'e_n');
M_.exo_names = char(M_.exo_names, 'e_p');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_p');
M_.exo_names_long = char(M_.exo_names_long, 'e_p');
M_.exo_names = char(M_.exo_names, 'e_w');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_w');
M_.exo_names_long = char(M_.exo_names_long, 'e_w');
M_.exo_names = char(M_.exo_names, 'e_e');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_e');
M_.exo_names_long = char(M_.exo_names_long, 'e_e');
M_.exo_names = char(M_.exo_names, 'e_r');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_r');
M_.exo_names_long = char(M_.exo_names_long, 'e_r');
M_.endo_names = 'uh';
M_.endo_names_tex = 'uh';
M_.endo_names_long = 'uh';
M_.endo_names = char(M_.endo_names, 'uc');
M_.endo_names_tex = char(M_.endo_names_tex, 'uc');
M_.endo_names_long = char(M_.endo_names_long, 'uc');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'rr');
M_.endo_names_tex = char(M_.endo_names_tex, 'rr');
M_.endo_names_long = char(M_.endo_names_long, 'rr');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'v');
M_.endo_names_tex = char(M_.endo_names_tex, 'v');
M_.endo_names_long = char(M_.endo_names_long, 'v');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'de');
M_.endo_names_tex = char(M_.endo_names_tex, 'de');
M_.endo_names_long = char(M_.endo_names_long, 'de');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'ne');
M_.endo_names_tex = char(M_.endo_names_tex, 'ne');
M_.endo_names_long = char(M_.endo_names_long, 'ne');
M_.endo_names = char(M_.endo_names, 'mk');
M_.endo_names_tex = char(M_.endo_names_tex, 'mk');
M_.endo_names_long = char(M_.endo_names_long, 'mk');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'Psi');
M_.endo_names_tex = char(M_.endo_names_tex, 'Psi');
M_.endo_names_long = char(M_.endo_names_long, 'Psi');
M_.endo_names = char(M_.endo_names, 'pic');
M_.endo_names_tex = char(M_.endo_names_tex, 'pic');
M_.endo_names_long = char(M_.endo_names_long, 'pic');
M_.endo_names = char(M_.endo_names, 'p');
M_.endo_names_tex = char(M_.endo_names_tex, 'p');
M_.endo_names_long = char(M_.endo_names_long, 'p');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'yd');
M_.endo_names_tex = char(M_.endo_names_tex, 'yd');
M_.endo_names_long = char(M_.endo_names_long, 'yd');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'he');
M_.endo_names_tex = char(M_.endo_names_tex, 'he');
M_.endo_names_long = char(M_.endo_names_long, 'he');
M_.endo_names = char(M_.endo_names, 'hc');
M_.endo_names_tex = char(M_.endo_names_tex, 'hc');
M_.endo_names_long = char(M_.endo_names_long, 'hc');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'u');
M_.endo_names_tex = char(M_.endo_names_tex, 'u');
M_.endo_names_long = char(M_.endo_names_long, 'u');
M_.endo_names = char(M_.endo_names, 'Phiu');
M_.endo_names_tex = char(M_.endo_names_tex, 'Phiu');
M_.endo_names_long = char(M_.endo_names_long, 'Phiu');
M_.endo_names = char(M_.endo_names, 'ku');
M_.endo_names_tex = char(M_.endo_names_tex, 'ku');
M_.endo_names_long = char(M_.endo_names_long, 'ku');
M_.endo_names = char(M_.endo_names, 'rK');
M_.endo_names_tex = char(M_.endo_names_tex, 'rK');
M_.endo_names_long = char(M_.endo_names_long, 'rK');
M_.endo_names = char(M_.endo_names, 'rL');
M_.endo_names_tex = char(M_.endo_names_tex, 'rL');
M_.endo_names_long = char(M_.endo_names_long, 'rL');
M_.endo_names = char(M_.endo_names, 'l');
M_.endo_names_tex = char(M_.endo_names_tex, 'l');
M_.endo_names_long = char(M_.endo_names_long, 'l');
M_.endo_names = char(M_.endo_names, 'nn');
M_.endo_names_tex = char(M_.endo_names_tex, 'nn');
M_.endo_names_long = char(M_.endo_names_long, 'nn');
M_.endo_names = char(M_.endo_names, 's');
M_.endo_names_tex = char(M_.endo_names_tex, 's');
M_.endo_names_long = char(M_.endo_names_long, 's');
M_.endo_names = char(M_.endo_names, 'omega');
M_.endo_names_tex = char(M_.endo_names_tex, 'omega');
M_.endo_names_long = char(M_.endo_names_long, 'omega');
M_.endo_names = char(M_.endo_names, 'mcL');
M_.endo_names_tex = char(M_.endo_names_tex, 'mcL');
M_.endo_names_long = char(M_.endo_names_long, 'mcL');
M_.endo_names = char(M_.endo_names, 'eta');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta');
M_.endo_names_long = char(M_.endo_names_long, 'eta');
M_.endo_names = char(M_.endo_names, 'w_sup');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_sup');
M_.endo_names_long = char(M_.endo_names_long, 'w_sup');
M_.endo_names = char(M_.endo_names, 'w_inf');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_inf');
M_.endo_names_long = char(M_.endo_names_long, 'w_inf');
M_.endo_names = char(M_.endo_names, 'mut_w');
M_.endo_names_tex = char(M_.endo_names_tex, 'mut\_w');
M_.endo_names_long = char(M_.endo_names_long, 'mut_w');
M_.endo_names = char(M_.endo_names, 'wh');
M_.endo_names_tex = char(M_.endo_names_tex, 'wh');
M_.endo_names_long = char(M_.endo_names_long, 'wh');
M_.endo_names = char(M_.endo_names, 'pi_w');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_w');
M_.endo_names_long = char(M_.endo_names_long, 'pi_w');
M_.endo_names = char(M_.endo_names, 'mut_p');
M_.endo_names_tex = char(M_.endo_names_tex, 'mut\_p');
M_.endo_names_long = char(M_.endo_names_long, 'mut_p');
M_.endo_names = char(M_.endo_names, 'mut_L');
M_.endo_names_tex = char(M_.endo_names_tex, 'mut\_L');
M_.endo_names_long = char(M_.endo_names_long, 'mut_L');
M_.endo_names = char(M_.endo_names, 's_a');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_a');
M_.endo_names_long = char(M_.endo_names_long, 's_a');
M_.endo_names = char(M_.endo_names, 's_g');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_g');
M_.endo_names_long = char(M_.endo_names_long, 's_g');
M_.endo_names = char(M_.endo_names, 's_b');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_b');
M_.endo_names_long = char(M_.endo_names_long, 's_b');
M_.endo_names = char(M_.endo_names, 's_i');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_i');
M_.endo_names_long = char(M_.endo_names_long, 's_i');
M_.endo_names = char(M_.endo_names, 's_l');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_l');
M_.endo_names_long = char(M_.endo_names_long, 's_l');
M_.endo_names = char(M_.endo_names, 's_n');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_n');
M_.endo_names_long = char(M_.endo_names_long, 's_n');
M_.endo_names = char(M_.endo_names, 's_p');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_p');
M_.endo_names_long = char(M_.endo_names_long, 's_p');
M_.endo_names = char(M_.endo_names, 's_w');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_w');
M_.endo_names_long = char(M_.endo_names_long, 's_w');
M_.endo_names = char(M_.endo_names, 's_e');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_e');
M_.endo_names_long = char(M_.endo_names_long, 's_e');
M_.endo_names = char(M_.endo_names, 's_r');
M_.endo_names_tex = char(M_.endo_names_tex, 's\_r');
M_.endo_names_long = char(M_.endo_names_long, 's_r');
M_.endo_names = char(M_.endo_names, 'ln_yd');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_yd');
M_.endo_names_long = char(M_.endo_names_long, 'ln_yd');
M_.endo_names = char(M_.endo_names, 'ln_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_y');
M_.endo_names_long = char(M_.endo_names_long, 'ln_y');
M_.endo_names = char(M_.endo_names, 'ln_c');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_c');
M_.endo_names_long = char(M_.endo_names_long, 'ln_c');
M_.endo_names = char(M_.endo_names, 'ln_n');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_n');
M_.endo_names_long = char(M_.endo_names_long, 'ln_n');
M_.endo_names = char(M_.endo_names, 'ln_r');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_r');
M_.endo_names_long = char(M_.endo_names_long, 'ln_r');
M_.endo_names = char(M_.endo_names, 'ln_pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_pi');
M_.endo_names_long = char(M_.endo_names_long, 'ln_pi');
M_.endo_names = char(M_.endo_names, 'ln_i');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_i');
M_.endo_names_long = char(M_.endo_names_long, 'ln_i');
M_.endo_names = char(M_.endo_names, 'ln_ne');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_ne');
M_.endo_names_long = char(M_.endo_names_long, 'ln_ne');
M_.endo_names = char(M_.endo_names, 'ln_rL');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_rL');
M_.endo_names_long = char(M_.endo_names_long, 'ln_rL');
M_.endo_names = char(M_.endo_names, 'ln_l');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_l');
M_.endo_names_long = char(M_.endo_names_long, 'ln_l');
M_.endo_names = char(M_.endo_names, 'ln_p');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_p');
M_.endo_names_long = char(M_.endo_names_long, 'ln_p');
M_.endo_names = char(M_.endo_names, 'ln_v');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_v');
M_.endo_names_long = char(M_.endo_names_long, 'ln_v');
M_.endo_names = char(M_.endo_names, 'y_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'y_obs');
M_.endo_names = char(M_.endo_names, 'c_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'c\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'c_obs');
M_.endo_names = char(M_.endo_names, 'i_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'i\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'i_obs');
M_.endo_names = char(M_.endo_names, 'pi_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'pi_obs');
M_.endo_names = char(M_.endo_names, 'r_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'r_obs');
M_.endo_names = char(M_.endo_names, 'w_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'w_obs');
M_.endo_names = char(M_.endo_names, 'h_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'h\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'h_obs');
M_.endo_names = char(M_.endo_names, 'ne_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'ne\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'ne_obs');
M_.endo_names = char(M_.endo_names, 'rL_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'rL\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'rL_obs');
M_.endo_names = char(M_.endo_names, 'l_obs');
M_.endo_names_tex = char(M_.endo_names_tex, 'l\_obs');
M_.endo_names_long = char(M_.endo_names_long, 'l_obs');
M_.endo_names = char(M_.endo_names, 'uhn');
M_.endo_names_tex = char(M_.endo_names_tex, 'uhn');
M_.endo_names_long = char(M_.endo_names_long, 'uhn');
M_.endo_names = char(M_.endo_names, 'ucn');
M_.endo_names_tex = char(M_.endo_names_tex, 'ucn');
M_.endo_names_long = char(M_.endo_names_long, 'ucn');
M_.endo_names = char(M_.endo_names, 'wn');
M_.endo_names_tex = char(M_.endo_names_tex, 'wn');
M_.endo_names_long = char(M_.endo_names_long, 'wn');
M_.endo_names = char(M_.endo_names, 'rrn');
M_.endo_names_tex = char(M_.endo_names_tex, 'rrn');
M_.endo_names_long = char(M_.endo_names_long, 'rrn');
M_.endo_names = char(M_.endo_names, 'vn');
M_.endo_names_tex = char(M_.endo_names_tex, 'vn');
M_.endo_names_long = char(M_.endo_names_long, 'vn');
M_.endo_names = char(M_.endo_names, 'dn');
M_.endo_names_tex = char(M_.endo_names_tex, 'dn');
M_.endo_names_long = char(M_.endo_names_long, 'dn');
M_.endo_names = char(M_.endo_names, 'den');
M_.endo_names_tex = char(M_.endo_names_tex, 'den');
M_.endo_names_long = char(M_.endo_names_long, 'den');
M_.endo_names = char(M_.endo_names, 'nnn');
M_.endo_names_tex = char(M_.endo_names_tex, 'nnn');
M_.endo_names_long = char(M_.endo_names_long, 'nnn');
M_.endo_names = char(M_.endo_names, 'nen');
M_.endo_names_tex = char(M_.endo_names_tex, 'nen');
M_.endo_names_long = char(M_.endo_names_long, 'nen');
M_.endo_names = char(M_.endo_names, 'mcn');
M_.endo_names_tex = char(M_.endo_names_tex, 'mcn');
M_.endo_names_long = char(M_.endo_names_long, 'mcn');
M_.endo_names = char(M_.endo_names, 'pn');
M_.endo_names_tex = char(M_.endo_names_tex, 'pn');
M_.endo_names_long = char(M_.endo_names_long, 'pn');
M_.endo_names = char(M_.endo_names, 'cn');
M_.endo_names_tex = char(M_.endo_names_tex, 'cn');
M_.endo_names_long = char(M_.endo_names_long, 'cn');
M_.endo_names = char(M_.endo_names, 'yn');
M_.endo_names_tex = char(M_.endo_names_tex, 'yn');
M_.endo_names_long = char(M_.endo_names_long, 'yn');
M_.endo_names = char(M_.endo_names, 'ydn');
M_.endo_names_tex = char(M_.endo_names_tex, 'ydn');
M_.endo_names_long = char(M_.endo_names_long, 'ydn');
M_.endo_names = char(M_.endo_names, 'hn');
M_.endo_names_tex = char(M_.endo_names_tex, 'hn');
M_.endo_names_long = char(M_.endo_names_long, 'hn');
M_.endo_names = char(M_.endo_names, 'hen');
M_.endo_names_tex = char(M_.endo_names_tex, 'hen');
M_.endo_names_long = char(M_.endo_names_long, 'hen');
M_.endo_names = char(M_.endo_names, 'hcn');
M_.endo_names_tex = char(M_.endo_names_tex, 'hcn');
M_.endo_names_long = char(M_.endo_names_long, 'hcn');
M_.endo_names = char(M_.endo_names, 'zn');
M_.endo_names_tex = char(M_.endo_names_tex, 'zn');
M_.endo_names_long = char(M_.endo_names_long, 'zn');
M_.endo_names = char(M_.endo_names, 'kn');
M_.endo_names_tex = char(M_.endo_names_tex, 'kn');
M_.endo_names_long = char(M_.endo_names_long, 'kn');
M_.endo_names = char(M_.endo_names, 'in');
M_.endo_names_tex = char(M_.endo_names_tex, 'in');
M_.endo_names_long = char(M_.endo_names_long, 'in');
M_.endo_names = char(M_.endo_names, 'qn');
M_.endo_names_tex = char(M_.endo_names_tex, 'qn');
M_.endo_names_long = char(M_.endo_names_long, 'qn');
M_.endo_names = char(M_.endo_names, 'un');
M_.endo_names_tex = char(M_.endo_names_tex, 'un');
M_.endo_names_long = char(M_.endo_names_long, 'un');
M_.endo_names = char(M_.endo_names, 'Phiun');
M_.endo_names_tex = char(M_.endo_names_tex, 'Phiun');
M_.endo_names_long = char(M_.endo_names_long, 'Phiun');
M_.endo_names = char(M_.endo_names, 'kun');
M_.endo_names_tex = char(M_.endo_names_tex, 'kun');
M_.endo_names_long = char(M_.endo_names_long, 'kun');
M_.endo_names = char(M_.endo_names, 'rKn');
M_.endo_names_tex = char(M_.endo_names_tex, 'rKn');
M_.endo_names_long = char(M_.endo_names_long, 'rKn');
M_.endo_names = char(M_.endo_names, 'rLn');
M_.endo_names_tex = char(M_.endo_names_tex, 'rLn');
M_.endo_names_long = char(M_.endo_names_long, 'rLn');
M_.endo_names = char(M_.endo_names, 'lnn');
M_.endo_names_tex = char(M_.endo_names_tex, 'lnn');
M_.endo_names_long = char(M_.endo_names_long, 'lnn');
M_.endo_names = char(M_.endo_names, 'nnnn');
M_.endo_names_tex = char(M_.endo_names_tex, 'nnnn');
M_.endo_names_long = char(M_.endo_names_long, 'nnnn');
M_.endo_names = char(M_.endo_names, 'sn');
M_.endo_names_tex = char(M_.endo_names_tex, 'sn');
M_.endo_names_long = char(M_.endo_names_long, 'sn');
M_.endo_names = char(M_.endo_names, 'omegan');
M_.endo_names_tex = char(M_.endo_names_tex, 'omegan');
M_.endo_names_long = char(M_.endo_names_long, 'omegan');
M_.endo_names = char(M_.endo_names, 'mcLn');
M_.endo_names_tex = char(M_.endo_names_tex, 'mcLn');
M_.endo_names_long = char(M_.endo_names_long, 'mcLn');
M_.endo_names = char(M_.endo_names, 'etan');
M_.endo_names_tex = char(M_.endo_names_tex, 'etan');
M_.endo_names_long = char(M_.endo_names_long, 'etan');
M_.endo_names = char(M_.endo_names, 'w_supn');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_supn');
M_.endo_names_long = char(M_.endo_names_long, 'w_supn');
M_.endo_names = char(M_.endo_names, 'w_infn');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_infn');
M_.endo_names_long = char(M_.endo_names_long, 'w_infn');
M_.endo_names = char(M_.endo_names, 'mut_wn');
M_.endo_names_tex = char(M_.endo_names_tex, 'mut\_wn');
M_.endo_names_long = char(M_.endo_names_long, 'mut_wn');
M_.endo_names = char(M_.endo_names, 'whn');
M_.endo_names_tex = char(M_.endo_names_tex, 'whn');
M_.endo_names_long = char(M_.endo_names_long, 'whn');
M_.endo_names = char(M_.endo_names, 'pi_wn');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_wn');
M_.endo_names_long = char(M_.endo_names_long, 'pi_wn');
M_.endo_names = char(M_.endo_names, 'mut_Ln');
M_.endo_names_tex = char(M_.endo_names_tex, 'mut\_Ln');
M_.endo_names_long = char(M_.endo_names_long, 'mut_Ln');
M_.endo_names = char(M_.endo_names, 'ln_yn');
M_.endo_names_tex = char(M_.endo_names_tex, 'ln\_yn');
M_.endo_names_long = char(M_.endo_names_long, 'ln_yn');
M_.endo_names = char(M_.endo_names, 'AUX_EXO_LAG_121_0');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_EXO\_LAG\_121\_0');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_EXO_LAG_121_0');
M_.endo_names = char(M_.endo_names, 'AUX_EXO_LAG_122_0');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_EXO\_LAG\_122\_0');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_EXO_LAG_122_0');
M_.endo_partitions = struct();
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'mu');
M_.param_names_tex = char(M_.param_names_tex, 'mu');
M_.param_names_long = char(M_.param_names_long, 'mu');
M_.param_names = char(M_.param_names, 'sigmaL');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaL');
M_.param_names_long = char(M_.param_names_long, 'sigmaL');
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'Fe');
M_.param_names_tex = char(M_.param_names_tex, 'Fe');
M_.param_names_long = char(M_.param_names_long, 'Fe');
M_.param_names = char(M_.param_names, 'kappa_P');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_P');
M_.param_names_long = char(M_.param_names_long, 'kappa_P');
M_.param_names = char(M_.param_names, 'xi_P');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_P');
M_.param_names_long = char(M_.param_names_long, 'xi_P');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'gy');
M_.param_names_tex = char(M_.param_names_tex, 'gy');
M_.param_names_long = char(M_.param_names_long, 'gy');
M_.param_names = char(M_.param_names, 'hh');
M_.param_names_tex = char(M_.param_names_tex, 'hh');
M_.param_names_long = char(M_.param_names_long, 'hh');
M_.param_names = char(M_.param_names, 'chi_I');
M_.param_names_tex = char(M_.param_names_tex, 'chi\_I');
M_.param_names_long = char(M_.param_names_long, 'chi_I');
M_.param_names = char(M_.param_names, 'chi_E');
M_.param_names_tex = char(M_.param_names_tex, 'chi\_E');
M_.param_names_long = char(M_.param_names_long, 'chi_E');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pi');
M_.param_names_long = char(M_.param_names_long, 'phi_pi');
M_.param_names = char(M_.param_names, 'phi_pic');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pic');
M_.param_names_long = char(M_.param_names_long, 'phi_pic');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_y');
M_.param_names_long = char(M_.param_names_long, 'phi_y');
M_.param_names = char(M_.param_names, 'phi_dy');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_dy');
M_.param_names_long = char(M_.param_names_long, 'phi_dy');
M_.param_names = char(M_.param_names, 'L_QK');
M_.param_names_tex = char(M_.param_names_tex, 'L\_QK');
M_.param_names_long = char(M_.param_names_long, 'L_QK');
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names_long = char(M_.param_names_long, 'R');
M_.param_names = char(M_.param_names, 'eta_d');
M_.param_names_tex = char(M_.param_names_tex, 'eta\_d');
M_.param_names_long = char(M_.param_names_long, 'eta_d');
M_.param_names = char(M_.param_names, 'H');
M_.param_names_tex = char(M_.param_names_tex, 'H');
M_.param_names_long = char(M_.param_names_long, 'H');
M_.param_names = char(M_.param_names, 'RL');
M_.param_names_tex = char(M_.param_names_tex, 'RL');
M_.param_names_long = char(M_.param_names_long, 'RL');
M_.param_names = char(M_.param_names, 'sigmaC');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaC');
M_.param_names_long = char(M_.param_names_long, 'sigmaC');
M_.param_names = char(M_.param_names, 'obsFactor');
M_.param_names_tex = char(M_.param_names_tex, 'obsFactor');
M_.param_names_long = char(M_.param_names_long, 'obsFactor');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'wmin');
M_.param_names_tex = char(M_.param_names_tex, 'wmin');
M_.param_names_long = char(M_.param_names_long, 'wmin');
M_.param_names = char(M_.param_names, 'ka');
M_.param_names_tex = char(M_.param_names_tex, 'ka');
M_.param_names_long = char(M_.param_names_long, 'ka');
M_.param_names = char(M_.param_names, 'TT');
M_.param_names_tex = char(M_.param_names_tex, 'TT');
M_.param_names_long = char(M_.param_names_long, 'TT');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'N_K');
M_.param_names_tex = char(M_.param_names_tex, 'N\_K');
M_.param_names_long = char(M_.param_names_long, 'N_K');
M_.param_names = char(M_.param_names, 'varkappa');
M_.param_names_tex = char(M_.param_names_tex, 'varkappa');
M_.param_names_long = char(M_.param_names_long, 'varkappa');
M_.param_names = char(M_.param_names, 'kappa_L');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_L');
M_.param_names_long = char(M_.param_names_long, 'kappa_L');
M_.param_names = char(M_.param_names, 'mu_B');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_B');
M_.param_names_long = char(M_.param_names_long, 'mu_B');
M_.param_names = char(M_.param_names, 'mu_L');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_L');
M_.param_names_long = char(M_.param_names_long, 'mu_L');
M_.param_names = char(M_.param_names, 'mu_P');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_P');
M_.param_names_long = char(M_.param_names_long, 'mu_P');
M_.param_names = char(M_.param_names, 'mu_W');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_W');
M_.param_names_long = char(M_.param_names_long, 'mu_W');
M_.param_names = char(M_.param_names, 'kappa_W');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_W');
M_.param_names_long = char(M_.param_names_long, 'kappa_W');
M_.param_names = char(M_.param_names, 'xi_W');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_W');
M_.param_names_long = char(M_.param_names_long, 'xi_W');
M_.param_names = char(M_.param_names, 'u_p');
M_.param_names_tex = char(M_.param_names_tex, 'u\_p');
M_.param_names_long = char(M_.param_names_long, 'u_p');
M_.param_names = char(M_.param_names, 'u_w');
M_.param_names_tex = char(M_.param_names_tex, 'u\_w');
M_.param_names_long = char(M_.param_names_long, 'u_w');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'rho_g');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_g');
M_.param_names_long = char(M_.param_names_long, 'rho_g');
M_.param_names = char(M_.param_names, 'rho_b');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_b');
M_.param_names_long = char(M_.param_names_long, 'rho_b');
M_.param_names = char(M_.param_names, 'rho_i');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_i');
M_.param_names_long = char(M_.param_names_long, 'rho_i');
M_.param_names = char(M_.param_names, 'rho_l');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_l');
M_.param_names_long = char(M_.param_names_long, 'rho_l');
M_.param_names = char(M_.param_names, 'rho_n');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_n');
M_.param_names_long = char(M_.param_names_long, 'rho_n');
M_.param_names = char(M_.param_names, 'rho_p');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_p');
M_.param_names_long = char(M_.param_names_long, 'rho_p');
M_.param_names = char(M_.param_names, 'rho_w');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_w');
M_.param_names_long = char(M_.param_names_long, 'rho_w');
M_.param_names = char(M_.param_names, 'rho_e');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_e');
M_.param_names_long = char(M_.param_names_long, 'rho_e');
M_.param_names = char(M_.param_names, 'rho_r');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_r');
M_.param_names_long = char(M_.param_names_long, 'rho_r');
M_.param_names = char(M_.param_names, 'rho_ag');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_ag');
M_.param_names_long = char(M_.param_names_long, 'rho_ag');
M_.param_names = char(M_.param_names, 'mkn');
M_.param_names_tex = char(M_.param_names_tex, 'mkn');
M_.param_names_long = char(M_.param_names_long, 'mkn');
M_.param_names = char(M_.param_names, 'Psin');
M_.param_names_tex = char(M_.param_names_tex, 'Psin');
M_.param_names_long = char(M_.param_names_long, 'Psin');
M_.param_names = char(M_.param_names, 'rn');
M_.param_names_tex = char(M_.param_names_tex, 'rn');
M_.param_names_long = char(M_.param_names_long, 'rn');
M_.param_names = char(M_.param_names, 'pin');
M_.param_names_tex = char(M_.param_names_tex, 'pin');
M_.param_names_long = char(M_.param_names_long, 'pin');
M_.param_names = char(M_.param_names, 'picn');
M_.param_names_tex = char(M_.param_names_tex, 'picn');
M_.param_names_long = char(M_.param_names_long, 'picn');
M_.param_names = char(M_.param_names, 'mut_pn');
M_.param_names_tex = char(M_.param_names_tex, 'mut\_pn');
M_.param_names_long = char(M_.param_names_long, 'mut_pn');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 10;
M_.endo_nbr = 117;
M_.param_nbr = 61;
M_.orig_endo_nbr = 115;
M_.aux_vars(1).endo_index = 116;
M_.aux_vars(1).type = 3;
M_.aux_vars(1).orig_index = 7;
M_.aux_vars(1).orig_lead_lag = 0;
M_.aux_vars(2).endo_index = 117;
M_.aux_vars(2).type = 3;
M_.aux_vars(2).orig_index = 8;
M_.aux_vars(2).orig_lead_lag = 0;
M_.Sigma_e = zeros(10, 10);
M_.Correlation_matrix = eye(10, 10);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('US_PV15_rep_static');
erase_compiled_function('US_PV15_rep_dynamic');
M_.orig_eq_nbr = 115;
M_.eq_nbr = 117;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 39 0;
 0 40 156;
 1 41 157;
 2 42 0;
 0 43 0;
 3 44 158;
 0 45 159;
 0 46 160;
 0 47 161;
 4 48 0;
 5 49 162;
 0 50 0;
 0 51 0;
 0 52 0;
 6 53 163;
 7 54 0;
 8 55 0;
 9 56 164;
 10 57 0;
 0 58 0;
 0 59 0;
 0 60 0;
 0 61 0;
 11 62 0;
 12 63 165;
 13 64 166;
 0 65 0;
 0 66 0;
 0 67 0;
 0 68 167;
 14 69 168;
 15 70 169;
 0 71 0;
 16 72 0;
 0 73 170;
 0 74 0;
 0 75 171;
 0 76 0;
 0 77 172;
 0 78 0;
 0 79 0;
 0 80 173;
 0 81 0;
 0 82 0;
 17 83 0;
 18 84 0;
 19 85 0;
 20 86 174;
 21 87 0;
 22 88 0;
 23 89 0;
 24 90 0;
 25 91 175;
 26 92 0;
 0 93 0;
 0 94 0;
 0 95 0;
 0 96 0;
 0 97 0;
 0 98 0;
 0 99 0;
 0 100 0;
 0 101 0;
 0 102 0;
 0 103 0;
 0 104 0;
 0 105 0;
 0 106 0;
 0 107 0;
 0 108 0;
 0 109 0;
 0 110 0;
 0 111 0;
 0 112 0;
 0 113 0;
 0 114 0;
 0 115 0;
 0 116 176;
 27 117 0;
 0 118 0;
 0 119 177;
 0 120 178;
 0 121 179;
 28 122 0;
 29 123 180;
 0 124 0;
 0 125 0;
 30 126 0;
 0 127 0;
 0 128 0;
 0 129 0;
 0 130 0;
 0 131 0;
 0 132 0;
 31 133 0;
 32 134 181;
 33 135 182;
 0 136 0;
 0 137 0;
 0 138 0;
 0 139 183;
 34 140 0;
 35 141 0;
 0 142 0;
 36 143 0;
 0 144 184;
 0 145 0;
 0 146 185;
 0 147 0;
 0 148 186;
 0 149 0;
 0 150 0;
 0 151 0;
 0 152 0;
 0 153 0;
 37 154 0;
 38 155 0;]';
M_.nstatic = 62;
M_.nfwrd   = 17;
M_.npred   = 24;
M_.nboth   = 14;
M_.nsfwrd   = 31;
M_.nspred   = 38;
M_.ndynamic   = 55;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:10];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(117, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(10, 1);
M_.params = NaN(61, 1);
M_.NNZDerivatives = [452; -1; -1];
close all;
M_.params( 1 ) = 0.992;
beta = M_.params( 1 );
M_.params( 3 ) = 0.025;
delta = M_.params( 3 );
M_.params( 10 ) = 3.8;
epsilon = M_.params( 10 );
M_.params( 4 ) = M_.params(10)/(M_.params(10)-1);
mu = M_.params( 4 );
M_.params( 7 ) = 5;
Fe = M_.params( 7 );
M_.params( 2 ) = .4;
alpha = M_.params( 2 );
M_.params( 11 ) = 0.18;
gy = M_.params( 11 );
M_.params( 40 ) = 1.5;
mu_W = M_.params( 40 );
M_.params( 34 ) = 0.5;
N_K = M_.params( 34 );
M_.params( 21 ) = 0.5;
L_QK = M_.params( 21 );
M_.params( 23 ) = 0.005;
eta_d = M_.params( 23 );
M_.params( 24 ) = 0.3333333333333333;
H = M_.params( 24 );
M_.params( 28 ) = 4;
theta = M_.params( 28 );
M_.params( 29 ) = 1;
gamma = M_.params( 29 );
M_.params( 22 ) = 1/M_.params(1);
R = M_.params( 22 );
M_.params( 25 ) = 1.0098*M_.params(22);
RL = M_.params( 25 );
M_.params( 27 ) = 100;
obsFactor = M_.params( 27 );
M_.params( 56 ) = M_.params(39);
mkn = M_.params( 56 );
M_.params( 57 ) = 0;
Psin = M_.params( 57 );
M_.params( 58 ) = M_.params(22)-1;
rn = M_.params( 58 );
M_.params( 59 ) = 1;
pin = M_.params( 59 );
M_.params( 60 ) = 1;
picn = M_.params( 60 );
M_.params( 61 ) = M_.params(39);
mut_pn = M_.params( 61 );
M_.params( 45 ) = 0.989438995688241;
rho_a = M_.params( 45 );
M_.params( 46 ) = 0.964501329345530;
rho_g = M_.params( 46 );
M_.params( 47 ) = 0.866065054296960;
rho_b = M_.params( 47 );
M_.params( 48 ) = 0.992049882657876;
rho_i = M_.params( 48 );
M_.params( 49 ) = 0.844755249127079;
rho_l = M_.params( 49 );
M_.params( 50 ) = 0.950156585064000;
rho_n = M_.params( 50 );
M_.params( 51 ) = 0.633888649873990;
rho_p = M_.params( 51 );
M_.params( 52 ) = 0.992861418199060;
rho_w = M_.params( 52 );
M_.params( 53 ) = 0.427852388654990;
rho_e = M_.params( 53 );
M_.params( 54 ) = 0.479821648248877;
rho_r = M_.params( 54 );
M_.params( 55 ) = 0.408305660498673;
rho_ag = M_.params( 55 );
M_.params( 26 ) = 1.412005386245682;
sigmaC = M_.params( 26 );
M_.params( 5 ) = 0.680598432322429;
sigmaL = M_.params( 5 );
M_.params( 41 ) = 57.815138488359928;
kappa_W = M_.params( 41 );
M_.params( 42 ) = 0.147527509676122;
xi_W = M_.params( 42 );
M_.params( 9 ) = 0.256415364089944;
xi_P = M_.params( 9 );
M_.params( 43 ) = 0.353493870232101;
u_p = M_.params( 43 );
M_.params( 44 ) = 0.880087358281597;
u_w = M_.params( 44 );
M_.params( 35 ) = 0.212721156496507;
varkappa = M_.params( 35 );
M_.params( 36 ) = 13.105140377417609;
kappa_L = M_.params( 36 );
M_.params( 37 ) = 0.074227073284695;
mu_B = M_.params( 37 );
M_.params( 8 ) = 52.900789993278337;
kappa_P = M_.params( 8 );
M_.params( 12 ) = 0.633600823023379;
hh = M_.params( 12 );
M_.params( 13 ) = 8.052744763552568;
chi_I = M_.params( 13 );
M_.params( 16 ) = 0.874082619958901;
rho = M_.params( 16 );
M_.params( 17 ) = 1.442040888486226;
phi_pi = M_.params( 17 );
M_.params( 18 ) = 0;
phi_pic = M_.params( 18 );
M_.params( 20 ) = 0.086084005670054;
phi_dy = M_.params( 20 );
M_.params( 19 ) = 0.107573728902499;
phi_y = M_.params( 19 );
M_.params( 15 ) = 0.865371316948408;
psi = M_.params( 15 );
M_.params( 14 ) = 0.911043517537755;
chi_E = M_.params( 14 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.96)^2;
M_.Sigma_e(2, 2) = (2.32)^2;
M_.Sigma_e(7, 7) = (0.05)^2;
M_.Sigma_e(10, 10) = (0.09)^2;
resid(1);
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 25;
options_.nograph = 1;
options_.order = 1;
var_list_ = char('ln_yd','ln_c','ln_i','ln_pi','ln_y','ln_n','ln_v','ln_r');
info = stoch_simul(var_list_);
save('US_PV15_rep_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('US_PV15_rep_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('US_PV15_rep_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('US_PV15_rep_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('US_PV15_rep_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('US_PV15_rep_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('US_PV15_rep_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
disp('Note: 2 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
