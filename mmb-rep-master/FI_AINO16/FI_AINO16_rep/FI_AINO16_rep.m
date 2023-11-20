%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'FI_AINO16_rep';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('FI_AINO16_rep.log');
M_.exo_names = 'epsMU';
M_.exo_names_tex = 'epsMU';
M_.exo_names_long = 'epsMU';
M_.exo_names = char(M_.exo_names, 'epsLAMBDAK');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMBDAK');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMBDAK');
M_.exo_names = char(M_.exo_names, 'epsLAMBDALT');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMBDALT');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMBDALT');
M_.exo_names = char(M_.exo_names, 'epsLAMBDACY');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMBDACY');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMBDACY');
M_.exo_names = char(M_.exo_names, 'epsLAMBDACM');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMBDACM');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMBDACM');
M_.exo_names = char(M_.exo_names, 'epsLAMBDAIY');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMBDAIY');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMBDAIY');
M_.exo_names = char(M_.exo_names, 'epsUPSILON');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsUPSILON');
M_.exo_names_long = char(M_.exo_names_long, 'epsUPSILON');
M_.exo_names = char(M_.exo_names, 'epsUPSILONMC');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsUPSILONMC');
M_.exo_names_long = char(M_.exo_names_long, 'epsUPSILONMC');
M_.exo_names = char(M_.exo_names, 'epsUPSILONX');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsUPSILONX');
M_.exo_names_long = char(M_.exo_names_long, 'epsUPSILONX');
M_.exo_names = char(M_.exo_names, 'epsZETACH');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsZETACH');
M_.exo_names_long = char(M_.exo_names_long, 'epsZETACH');
M_.exo_names = char(M_.exo_names, 'epsLAMW');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsLAMW');
M_.exo_names_long = char(M_.exo_names_long, 'epsLAMW');
M_.exo_names = char(M_.exo_names, 'epsZETAEUR');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsZETAEUR');
M_.exo_names_long = char(M_.exo_names_long, 'epsZETAEUR');
M_.exo_names = char(M_.exo_names, 'epsdS');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsdS');
M_.exo_names_long = char(M_.exo_names_long, 'epsdS');
M_.exo_names = char(M_.exo_names, 'epsIG');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsIG');
M_.exo_names_long = char(M_.exo_names_long, 'epsIG');
M_.exo_names = char(M_.exo_names, 'epshG');
M_.exo_names_tex = char(M_.exo_names_tex, 'epshG');
M_.exo_names_long = char(M_.exo_names_long, 'epshG');
M_.exo_names = char(M_.exo_names, 'epsGF');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsGF');
M_.exo_names_long = char(M_.exo_names_long, 'epsGF');
M_.exo_names = char(M_.exo_names, 'epsPOILS');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsPOILS');
M_.exo_names_long = char(M_.exo_names_long, 'epsPOILS');
M_.exo_names = char(M_.exo_names, 'epsPRAWS');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsPRAWS');
M_.exo_names_long = char(M_.exo_names_long, 'epsPRAWS');
M_.exo_names = char(M_.exo_names, 'epsMW');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsMW');
M_.exo_names_long = char(M_.exo_names_long, 'epsMW');
M_.exo_names = char(M_.exo_names, 'epsPIEW');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsPIEW');
M_.exo_names_long = char(M_.exo_names_long, 'epsPIEW');
M_.exo_names = char(M_.exo_names, 'epsXX');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsXX');
M_.exo_names_long = char(M_.exo_names_long, 'epsXX');
M_.exo_names = char(M_.exo_names, 'epsEPSB');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsEPSB');
M_.exo_names_long = char(M_.exo_names_long, 'epsEPSB');
M_.exo_names = char(M_.exo_names, 'epsnuB');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsnuB');
M_.exo_names_long = char(M_.exo_names_long, 'epsnuB');
M_.exo_names = char(M_.exo_names, 'epsBankCapital');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsBankCapital');
M_.exo_names_long = char(M_.exo_names_long, 'epsBankCapital');
M_.exo_names = char(M_.exo_names, 'epsrEUR');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsrEUR');
M_.exo_names_long = char(M_.exo_names_long, 'epsrEUR');
M_.endo_names = 'psi';
M_.endo_names_tex = 'psi';
M_.endo_names_long = 'psi';
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'iH');
M_.endo_names_tex = char(M_.endo_names_tex, 'iH');
M_.endo_names_long = char(M_.endo_names_long, 'iH');
M_.endo_names = char(M_.endo_names, 'cH');
M_.endo_names_tex = char(M_.endo_names_tex, 'cH');
M_.endo_names_long = char(M_.endo_names_long, 'cH');
M_.endo_names = char(M_.endo_names, 'rFI');
M_.endo_names_tex = char(M_.endo_names_tex, 'rFI');
M_.endo_names_long = char(M_.endo_names_long, 'rFI');
M_.endo_names = char(M_.endo_names, 'rEUR');
M_.endo_names_tex = char(M_.endo_names_tex, 'rEUR');
M_.endo_names_long = char(M_.endo_names_long, 'rEUR');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'wF');
M_.endo_names_tex = char(M_.endo_names_tex, 'wF');
M_.endo_names_long = char(M_.endo_names_long, 'wF');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'hF');
M_.endo_names_tex = char(M_.endo_names_tex, 'hF');
M_.endo_names_long = char(M_.endo_names_long, 'hF');
M_.endo_names = char(M_.endo_names, 'mcY');
M_.endo_names_tex = char(M_.endo_names_tex, 'mcY');
M_.endo_names_long = char(M_.endo_names_long, 'mcY');
M_.endo_names = char(M_.endo_names, 'pieY');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieY');
M_.endo_names_long = char(M_.endo_names_long, 'pieY');
M_.endo_names = char(M_.endo_names, 'ds');
M_.endo_names_tex = char(M_.endo_names_tex, 'ds');
M_.endo_names_long = char(M_.endo_names_long, 'ds');
M_.endo_names = char(M_.endo_names, 'yC');
M_.endo_names_tex = char(M_.endo_names_tex, 'yC');
M_.endo_names_long = char(M_.endo_names_long, 'yC');
M_.endo_names = char(M_.endo_names, 'mC');
M_.endo_names_tex = char(M_.endo_names_tex, 'mC');
M_.endo_names_long = char(M_.endo_names_long, 'mC');
M_.endo_names = char(M_.endo_names, 'pC');
M_.endo_names_tex = char(M_.endo_names_tex, 'pC');
M_.endo_names_long = char(M_.endo_names_long, 'pC');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'pO');
M_.endo_names_tex = char(M_.endo_names_tex, 'pO');
M_.endo_names_long = char(M_.endo_names_long, 'pO');
M_.endo_names = char(M_.endo_names, 'iT');
M_.endo_names_tex = char(M_.endo_names_tex, 'iT');
M_.endo_names_long = char(M_.endo_names_long, 'iT');
M_.endo_names = char(M_.endo_names, 'yI');
M_.endo_names_tex = char(M_.endo_names_tex, 'yI');
M_.endo_names_long = char(M_.endo_names_long, 'yI');
M_.endo_names = char(M_.endo_names, 'mI');
M_.endo_names_tex = char(M_.endo_names_tex, 'mI');
M_.endo_names_long = char(M_.endo_names_long, 'mI');
M_.endo_names = char(M_.endo_names, 'pI');
M_.endo_names_tex = char(M_.endo_names_tex, 'pI');
M_.endo_names_long = char(M_.endo_names_long, 'pI');
M_.endo_names = char(M_.endo_names, 'rK');
M_.endo_names_tex = char(M_.endo_names_tex, 'rK');
M_.endo_names_long = char(M_.endo_names_long, 'rK');
M_.endo_names = char(M_.endo_names, 'yX');
M_.endo_names_tex = char(M_.endo_names_tex, 'yX');
M_.endo_names_long = char(M_.endo_names_long, 'yX');
M_.endo_names = char(M_.endo_names, 'mX');
M_.endo_names_tex = char(M_.endo_names_tex, 'mX');
M_.endo_names_long = char(M_.endo_names_long, 'mX');
M_.endo_names = char(M_.endo_names, 'mcX');
M_.endo_names_tex = char(M_.endo_names_tex, 'mcX');
M_.endo_names_long = char(M_.endo_names_long, 'mcX');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'pieX');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieX');
M_.endo_names_long = char(M_.endo_names_long, 'pieX');
M_.endo_names = char(M_.endo_names, 'pX');
M_.endo_names_tex = char(M_.endo_names_tex, 'pX');
M_.endo_names_long = char(M_.endo_names_long, 'pX');
M_.endo_names = char(M_.endo_names, 'pieMC');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieMC');
M_.endo_names_long = char(M_.endo_names_long, 'pieMC');
M_.endo_names = char(M_.endo_names, 'pieM');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieM');
M_.endo_names_long = char(M_.endo_names_long, 'pieM');
M_.endo_names = char(M_.endo_names, 'tbY');
M_.endo_names_tex = char(M_.endo_names_tex, 'tbY');
M_.endo_names_long = char(M_.endo_names_long, 'tbY');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names_long = char(M_.endo_names_long, 'm');
M_.endo_names = char(M_.endo_names, 'pM');
M_.endo_names_tex = char(M_.endo_names_tex, 'pM');
M_.endo_names_long = char(M_.endo_names_long, 'pM');
M_.endo_names = char(M_.endo_names, 'ToT');
M_.endo_names_tex = char(M_.endo_names_tex, 'ToT');
M_.endo_names_long = char(M_.endo_names_long, 'ToT');
M_.endo_names = char(M_.endo_names, 'bstar');
M_.endo_names_tex = char(M_.endo_names_tex, 'bstar');
M_.endo_names_long = char(M_.endo_names_long, 'bstar');
M_.endo_names = char(M_.endo_names, 'astar');
M_.endo_names_tex = char(M_.endo_names_tex, 'astar');
M_.endo_names_long = char(M_.endo_names_long, 'astar');
M_.endo_names = char(M_.endo_names, 'pieC');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieC');
M_.endo_names_long = char(M_.endo_names_long, 'pieC');
M_.endo_names = char(M_.endo_names, 'pieI');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieI');
M_.endo_names_long = char(M_.endo_names_long, 'pieI');
M_.endo_names = char(M_.endo_names, 'rs');
M_.endo_names_tex = char(M_.endo_names_tex, 'rs');
M_.endo_names_long = char(M_.endo_names_long, 'rs');
M_.endo_names = char(M_.endo_names, 'trY');
M_.endo_names_tex = char(M_.endo_names_tex, 'trY');
M_.endo_names_long = char(M_.endo_names_long, 'trY');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'mu');
M_.endo_names_tex = char(M_.endo_names_tex, 'mu');
M_.endo_names_long = char(M_.endo_names_long, 'mu');
M_.endo_names = char(M_.endo_names, 'lamK');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamK');
M_.endo_names_long = char(M_.endo_names_long, 'lamK');
M_.endo_names = char(M_.endo_names, 'lamLT');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamLT');
M_.endo_names_long = char(M_.endo_names_long, 'lamLT');
M_.endo_names = char(M_.endo_names, 'lamCY');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamCY');
M_.endo_names_long = char(M_.endo_names_long, 'lamCY');
M_.endo_names = char(M_.endo_names, 'lamCM');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamCM');
M_.endo_names_long = char(M_.endo_names_long, 'lamCM');
M_.endo_names = char(M_.endo_names, 'lamIY');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamIY');
M_.endo_names_long = char(M_.endo_names_long, 'lamIY');
M_.endo_names = char(M_.endo_names, 'upsilon');
M_.endo_names_tex = char(M_.endo_names_tex, 'upsilon');
M_.endo_names_long = char(M_.endo_names_long, 'upsilon');
M_.endo_names = char(M_.endo_names, 'upsilonMC');
M_.endo_names_tex = char(M_.endo_names_tex, 'upsilonMC');
M_.endo_names_long = char(M_.endo_names_long, 'upsilonMC');
M_.endo_names = char(M_.endo_names, 'upsilonX');
M_.endo_names_tex = char(M_.endo_names_tex, 'upsilonX');
M_.endo_names_long = char(M_.endo_names_long, 'upsilonX');
M_.endo_names = char(M_.endo_names, 'zetaCH');
M_.endo_names_tex = char(M_.endo_names_tex, 'zetaCH');
M_.endo_names_long = char(M_.endo_names_long, 'zetaCH');
M_.endo_names = char(M_.endo_names, 'zetaEUR');
M_.endo_names_tex = char(M_.endo_names_tex, 'zetaEUR');
M_.endo_names_long = char(M_.endo_names_long, 'zetaEUR');
M_.endo_names = char(M_.endo_names, 'lamW');
M_.endo_names_tex = char(M_.endo_names_tex, 'lamW');
M_.endo_names_long = char(M_.endo_names_long, 'lamW');
M_.endo_names = char(M_.endo_names, 'iG');
M_.endo_names_tex = char(M_.endo_names_tex, 'iG');
M_.endo_names_long = char(M_.endo_names_long, 'iG');
M_.endo_names = char(M_.endo_names, 'hG');
M_.endo_names_tex = char(M_.endo_names_tex, 'hG');
M_.endo_names_long = char(M_.endo_names_long, 'hG');
M_.endo_names = char(M_.endo_names, 'cGF');
M_.endo_names_tex = char(M_.endo_names_tex, 'cGF');
M_.endo_names_long = char(M_.endo_names_long, 'cGF');
M_.endo_names = char(M_.endo_names, 'wG');
M_.endo_names_tex = char(M_.endo_names_tex, 'wG');
M_.endo_names_long = char(M_.endo_names_long, 'wG');
M_.endo_names = char(M_.endo_names, 'pOILS');
M_.endo_names_tex = char(M_.endo_names_tex, 'pOILS');
M_.endo_names_long = char(M_.endo_names_long, 'pOILS');
M_.endo_names = char(M_.endo_names, 'mW');
M_.endo_names_tex = char(M_.endo_names_tex, 'mW');
M_.endo_names_long = char(M_.endo_names_long, 'mW');
M_.endo_names = char(M_.endo_names, 'pieW');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieW');
M_.endo_names_long = char(M_.endo_names_long, 'pieW');
M_.endo_names = char(M_.endo_names, 'pieOILS');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieOILS');
M_.endo_names_long = char(M_.endo_names_long, 'pieOILS');
M_.endo_names = char(M_.endo_names, 'epsX');
M_.endo_names_tex = char(M_.endo_names_tex, 'epsX');
M_.endo_names_long = char(M_.endo_names_long, 'epsX');
M_.endo_names = char(M_.endo_names, 'pRAWS');
M_.endo_names_tex = char(M_.endo_names_tex, 'pRAWS');
M_.endo_names_long = char(M_.endo_names_long, 'pRAWS');
M_.endo_names = char(M_.endo_names, 'pieRAWS');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieRAWS');
M_.endo_names_long = char(M_.endo_names_long, 'pieRAWS');
M_.endo_names = char(M_.endo_names, 'rb');
M_.endo_names_tex = char(M_.endo_names_tex, 'rb');
M_.endo_names_long = char(M_.endo_names_long, 'rb');
M_.endo_names = char(M_.endo_names, 'nwe');
M_.endo_names_tex = char(M_.endo_names_tex, 'nwe');
M_.endo_names_long = char(M_.endo_names_long, 'nwe');
M_.endo_names = char(M_.endo_names, 'btot');
M_.endo_names_tex = char(M_.endo_names_tex, 'btot');
M_.endo_names_long = char(M_.endo_names_long, 'btot');
M_.endo_names = char(M_.endo_names, 'epsb');
M_.endo_names_tex = char(M_.endo_names_tex, 'epsb');
M_.endo_names_long = char(M_.endo_names_long, 'epsb');
M_.endo_names = char(M_.endo_names, 'lev_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'lev\_e');
M_.endo_names_long = char(M_.endo_names_long, 'lev_e');
M_.endo_names = char(M_.endo_names, 'by');
M_.endo_names_tex = char(M_.endo_names_tex, 'by');
M_.endo_names_long = char(M_.endo_names_long, 'by');
M_.endo_names = char(M_.endo_names, 'RB');
M_.endo_names_tex = char(M_.endo_names_tex, 'RB');
M_.endo_names_long = char(M_.endo_names_long, 'RB');
M_.endo_names = char(M_.endo_names, 'kbank');
M_.endo_names_tex = char(M_.endo_names_tex, 'kbank');
M_.endo_names_long = char(M_.endo_names_long, 'kbank');
M_.endo_names = char(M_.endo_names, 'nuB');
M_.endo_names_tex = char(M_.endo_names_tex, 'nuB');
M_.endo_names_long = char(M_.endo_names_long, 'nuB');
M_.endo_names = char(M_.endo_names, 'epsKB');
M_.endo_names_tex = char(M_.endo_names_tex, 'epsKB');
M_.endo_names_long = char(M_.endo_names_long, 'epsKB');
M_.endo_names = char(M_.endo_names, 'bankprofits');
M_.endo_names_tex = char(M_.endo_names_tex, 'bankprofits');
M_.endo_names_long = char(M_.endo_names_long, 'bankprofits');
M_.endo_names = char(M_.endo_names, 'deposits');
M_.endo_names_tex = char(M_.endo_names_tex, 'deposits');
M_.endo_names_long = char(M_.endo_names_long, 'deposits');
M_.endo_names = char(M_.endo_names, 'bka');
M_.endo_names_tex = char(M_.endo_names_tex, 'bka');
M_.endo_names_long = char(M_.endo_names_long, 'bka');
M_.endo_names = char(M_.endo_names, 'rwage');
M_.endo_names_tex = char(M_.endo_names_tex, 'rwage');
M_.endo_names_long = char(M_.endo_names_long, 'rwage');
M_.param_names = 'bC';
M_.param_names_tex = 'bC';
M_.param_names_long = 'bC';
M_.param_names = char(M_.param_names, 'ssMU');
M_.param_names_tex = char(M_.param_names_tex, 'ssMU');
M_.param_names_long = char(M_.param_names_long, 'ssMU');
M_.param_names = char(M_.param_names, 'ssTAXWR');
M_.param_names_tex = char(M_.param_names_tex, 'ssTAXWR');
M_.param_names_long = char(M_.param_names_long, 'ssTAXWR');
M_.param_names = char(M_.param_names, 'ssRS');
M_.param_names_tex = char(M_.param_names_tex, 'ssRS');
M_.param_names_long = char(M_.param_names_long, 'ssRS');
M_.param_names = char(M_.param_names, 'ssRPOILS');
M_.param_names_tex = char(M_.param_names_tex, 'ssRPOILS');
M_.param_names_long = char(M_.param_names_long, 'ssRPOILS');
M_.param_names = char(M_.param_names, 'ssTAXCR');
M_.param_names_tex = char(M_.param_names_tex, 'ssTAXCR');
M_.param_names_long = char(M_.param_names_long, 'ssTAXCR');
M_.param_names = char(M_.param_names, 'bet');
M_.param_names_tex = char(M_.param_names_tex, 'bet');
M_.param_names_long = char(M_.param_names_long, 'bet');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'lambdaW');
M_.param_names_tex = char(M_.param_names_tex, 'lambdaW');
M_.param_names_long = char(M_.param_names_long, 'lambdaW');
M_.param_names = char(M_.param_names, 'rhoY');
M_.param_names_tex = char(M_.param_names_tex, 'rhoY');
M_.param_names_long = char(M_.param_names_long, 'rhoY');
M_.param_names = char(M_.param_names, 'deltaY');
M_.param_names_tex = char(M_.param_names_tex, 'deltaY');
M_.param_names_long = char(M_.param_names_long, 'deltaY');
M_.param_names = char(M_.param_names, 'ssLAMBDAK');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDAK');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDAK');
M_.param_names = char(M_.param_names, 'ssLAMBDALT');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDALT');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDALT');
M_.param_names = char(M_.param_names, 'ssTAXKR');
M_.param_names_tex = char(M_.param_names_tex, 'ssTAXKR');
M_.param_names_long = char(M_.param_names_long, 'ssTAXKR');
M_.param_names = char(M_.param_names, 'gamI');
M_.param_names_tex = char(M_.param_names_tex, 'gamI');
M_.param_names_long = char(M_.param_names_long, 'gamI');
M_.param_names = char(M_.param_names, 'sspieY');
M_.param_names_tex = char(M_.param_names_tex, 'sspieY');
M_.param_names_long = char(M_.param_names_long, 'sspieY');
M_.param_names = char(M_.param_names, 'phia');
M_.param_names_tex = char(M_.param_names_tex, 'phia');
M_.param_names_long = char(M_.param_names_long, 'phia');
M_.param_names = char(M_.param_names, 'xiW');
M_.param_names_tex = char(M_.param_names_tex, 'xiW');
M_.param_names_long = char(M_.param_names_long, 'xiW');
M_.param_names = char(M_.param_names, 'sigmaL');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaL');
M_.param_names_long = char(M_.param_names_long, 'sigmaL');
M_.param_names = char(M_.param_names, 'ssHG');
M_.param_names_tex = char(M_.param_names_tex, 'ssHG');
M_.param_names_long = char(M_.param_names_long, 'ssHG');
M_.param_names = char(M_.param_names, 'ssTAXFR');
M_.param_names_tex = char(M_.param_names_tex, 'ssTAXFR');
M_.param_names_long = char(M_.param_names_long, 'ssTAXFR');
M_.param_names = char(M_.param_names, 'ssUPSILON');
M_.param_names_tex = char(M_.param_names_tex, 'ssUPSILON');
M_.param_names_long = char(M_.param_names_long, 'ssUPSILON');
M_.param_names = char(M_.param_names, 'ssUPSILONMC');
M_.param_names_tex = char(M_.param_names_tex, 'ssUPSILONMC');
M_.param_names_long = char(M_.param_names_long, 'ssUPSILONMC');
M_.param_names = char(M_.param_names, 'ssUPSILONMI');
M_.param_names_tex = char(M_.param_names_tex, 'ssUPSILONMI');
M_.param_names_long = char(M_.param_names_long, 'ssUPSILONMI');
M_.param_names = char(M_.param_names, 'ssUPSILONMX');
M_.param_names_tex = char(M_.param_names_tex, 'ssUPSILONMX');
M_.param_names_long = char(M_.param_names_long, 'ssUPSILONMX');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names_long = char(M_.param_names_long, 'zeta');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'deltaI');
M_.param_names_tex = char(M_.param_names_tex, 'deltaI');
M_.param_names_long = char(M_.param_names_long, 'deltaI');
M_.param_names = char(M_.param_names, 'deltaC');
M_.param_names_tex = char(M_.param_names_tex, 'deltaC');
M_.param_names_long = char(M_.param_names_long, 'deltaC');
M_.param_names = char(M_.param_names, 'ssLAMBDACY');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDACY');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDACY');
M_.param_names = char(M_.param_names, 'rhoC');
M_.param_names_tex = char(M_.param_names_tex, 'rhoC');
M_.param_names_long = char(M_.param_names_long, 'rhoC');
M_.param_names = char(M_.param_names, 'ssLAMBDACM');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDACM');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDACM');
M_.param_names = char(M_.param_names, 'ssGCF');
M_.param_names_tex = char(M_.param_names_tex, 'ssGCF');
M_.param_names_long = char(M_.param_names_long, 'ssGCF');
M_.param_names = char(M_.param_names, 'gamCM');
M_.param_names_tex = char(M_.param_names_tex, 'gamCM');
M_.param_names_long = char(M_.param_names_long, 'gamCM');
M_.param_names = char(M_.param_names, 'zetaMC');
M_.param_names_tex = char(M_.param_names_tex, 'zetaMC');
M_.param_names_long = char(M_.param_names_long, 'zetaMC');
M_.param_names = char(M_.param_names, 'ssLAMBDAIY');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDAIY');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDAIY');
M_.param_names = char(M_.param_names, 'ssLAMBDAIM');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDAIM');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDAIM');
M_.param_names = char(M_.param_names, 'rhoI');
M_.param_names_tex = char(M_.param_names_tex, 'rhoI');
M_.param_names_long = char(M_.param_names_long, 'rhoI');
M_.param_names = char(M_.param_names, 'gamIM');
M_.param_names_tex = char(M_.param_names_tex, 'gamIM');
M_.param_names_long = char(M_.param_names_long, 'gamIM');
M_.param_names = char(M_.param_names, 'ssIG');
M_.param_names_tex = char(M_.param_names_tex, 'ssIG');
M_.param_names_long = char(M_.param_names_long, 'ssIG');
M_.param_names = char(M_.param_names, 'rhoX');
M_.param_names_tex = char(M_.param_names_tex, 'rhoX');
M_.param_names_long = char(M_.param_names_long, 'rhoX');
M_.param_names = char(M_.param_names, 'deltaX');
M_.param_names_tex = char(M_.param_names_tex, 'deltaX');
M_.param_names_long = char(M_.param_names_long, 'deltaX');
M_.param_names = char(M_.param_names, 'thetaX');
M_.param_names_tex = char(M_.param_names_tex, 'thetaX');
M_.param_names_long = char(M_.param_names_long, 'thetaX');
M_.param_names = char(M_.param_names, 'zetaX');
M_.param_names_tex = char(M_.param_names_tex, 'zetaX');
M_.param_names_long = char(M_.param_names_long, 'zetaX');
M_.param_names = char(M_.param_names, 'ssLAMBDAXY');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDAXY');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDAXY');
M_.param_names = char(M_.param_names, 'ssLAMBDAXM');
M_.param_names_tex = char(M_.param_names_tex, 'ssLAMBDAXM');
M_.param_names_long = char(M_.param_names_long, 'ssLAMBDAXM');
M_.param_names = char(M_.param_names, 'ssUPSILONF');
M_.param_names_tex = char(M_.param_names_tex, 'ssUPSILONF');
M_.param_names_long = char(M_.param_names_long, 'ssUPSILONF');
M_.param_names = char(M_.param_names, 'thetaMC');
M_.param_names_tex = char(M_.param_names_tex, 'thetaMC');
M_.param_names_long = char(M_.param_names_long, 'thetaMC');
M_.param_names = char(M_.param_names, 'omegaMC');
M_.param_names_tex = char(M_.param_names_tex, 'omegaMC');
M_.param_names_long = char(M_.param_names_long, 'omegaMC');
M_.param_names = char(M_.param_names, 'rhoMU');
M_.param_names_tex = char(M_.param_names_tex, 'rhoMU');
M_.param_names_long = char(M_.param_names_long, 'rhoMU');
M_.param_names = char(M_.param_names, 'rhoLAMBDAK');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMBDAK');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMBDAK');
M_.param_names = char(M_.param_names, 'rhoLAMBDALT');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMBDALT');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMBDALT');
M_.param_names = char(M_.param_names, 'rhoLAMBDACY');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMBDACY');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMBDACY');
M_.param_names = char(M_.param_names, 'rhoLAMBDACM');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMBDACM');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMBDACM');
M_.param_names = char(M_.param_names, 'rhoLAMBDAIY');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMBDAIY');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMBDAIY');
M_.param_names = char(M_.param_names, 'rhoUPSILON');
M_.param_names_tex = char(M_.param_names_tex, 'rhoUPSILON');
M_.param_names_long = char(M_.param_names_long, 'rhoUPSILON');
M_.param_names = char(M_.param_names, 'rhoUPSILONMC');
M_.param_names_tex = char(M_.param_names_tex, 'rhoUPSILONMC');
M_.param_names_long = char(M_.param_names_long, 'rhoUPSILONMC');
M_.param_names = char(M_.param_names, 'rhoUPSILONX');
M_.param_names_tex = char(M_.param_names_tex, 'rhoUPSILONX');
M_.param_names_long = char(M_.param_names_long, 'rhoUPSILONX');
M_.param_names = char(M_.param_names, 'rhoZETACH');
M_.param_names_tex = char(M_.param_names_tex, 'rhoZETACH');
M_.param_names_long = char(M_.param_names_long, 'rhoZETACH');
M_.param_names = char(M_.param_names, 'rhoZETAEUR');
M_.param_names_tex = char(M_.param_names_tex, 'rhoZETAEUR');
M_.param_names_long = char(M_.param_names_long, 'rhoZETAEUR');
M_.param_names = char(M_.param_names, 'rhoLAMW');
M_.param_names_tex = char(M_.param_names_tex, 'rhoLAMW');
M_.param_names_long = char(M_.param_names_long, 'rhoLAMW');
M_.param_names = char(M_.param_names, 'rhoiG');
M_.param_names_tex = char(M_.param_names_tex, 'rhoiG');
M_.param_names_long = char(M_.param_names_long, 'rhoiG');
M_.param_names = char(M_.param_names, 'rhohG');
M_.param_names_tex = char(M_.param_names_tex, 'rhohG');
M_.param_names_long = char(M_.param_names_long, 'rhohG');
M_.param_names = char(M_.param_names, 'rhocFG');
M_.param_names_tex = char(M_.param_names_tex, 'rhocFG');
M_.param_names_long = char(M_.param_names_long, 'rhocFG');
M_.param_names = char(M_.param_names, 'rhoPOILS');
M_.param_names_tex = char(M_.param_names_tex, 'rhoPOILS');
M_.param_names_long = char(M_.param_names_long, 'rhoPOILS');
M_.param_names = char(M_.param_names, 'rhoPRAWS');
M_.param_names_tex = char(M_.param_names_tex, 'rhoPRAWS');
M_.param_names_long = char(M_.param_names_long, 'rhoPRAWS');
M_.param_names = char(M_.param_names, 'rhomW');
M_.param_names_tex = char(M_.param_names_tex, 'rhomW');
M_.param_names_long = char(M_.param_names_long, 'rhomW');
M_.param_names = char(M_.param_names, 'rhopieW');
M_.param_names_tex = char(M_.param_names_tex, 'rhopieW');
M_.param_names_long = char(M_.param_names_long, 'rhopieW');
M_.param_names = char(M_.param_names, 'rhodS');
M_.param_names_tex = char(M_.param_names_tex, 'rhodS');
M_.param_names_long = char(M_.param_names_long, 'rhodS');
M_.param_names = char(M_.param_names, 'rhoepsX');
M_.param_names_tex = char(M_.param_names_tex, 'rhoepsX');
M_.param_names_long = char(M_.param_names_long, 'rhoepsX');
M_.param_names = char(M_.param_names, 'omegaOIL');
M_.param_names_tex = char(M_.param_names_tex, 'omegaOIL');
M_.param_names_long = char(M_.param_names_long, 'omegaOIL');
M_.param_names = char(M_.param_names, 'omegaRAW');
M_.param_names_tex = char(M_.param_names_tex, 'omegaRAW');
M_.param_names_long = char(M_.param_names_long, 'omegaRAW');
M_.param_names = char(M_.param_names, 'rhorEUR');
M_.param_names_tex = char(M_.param_names_tex, 'rhorEUR');
M_.param_names_long = char(M_.param_names_long, 'rhorEUR');
M_.param_names = char(M_.param_names, 'RMCX');
M_.param_names_tex = char(M_.param_names_tex, 'RMCX');
M_.param_names_long = char(M_.param_names_long, 'RMCX');
M_.param_names = char(M_.param_names, 'sigmaW');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaW');
M_.param_names_long = char(M_.param_names_long, 'sigmaW');
M_.param_names = char(M_.param_names, 'WEHE');
M_.param_names_tex = char(M_.param_names_tex, 'WEHE');
M_.param_names_long = char(M_.param_names_long, 'WEHE');
M_.param_names = char(M_.param_names, 'ssnuBank');
M_.param_names_tex = char(M_.param_names_tex, 'ssnuBank');
M_.param_names_long = char(M_.param_names_long, 'ssnuBank');
M_.param_names = char(M_.param_names, 'kappaB');
M_.param_names_tex = char(M_.param_names_tex, 'kappaB');
M_.param_names_long = char(M_.param_names_long, 'kappaB');
M_.param_names = char(M_.param_names, 'kappaKB');
M_.param_names_tex = char(M_.param_names_tex, 'kappaKB');
M_.param_names_long = char(M_.param_names_long, 'kappaKB');
M_.param_names = char(M_.param_names, 'rhoepsB');
M_.param_names_tex = char(M_.param_names_tex, 'rhoepsB');
M_.param_names_long = char(M_.param_names_long, 'rhoepsB');
M_.param_names = char(M_.param_names, 'rhoepsKB');
M_.param_names_tex = char(M_.param_names_tex, 'rhoepsKB');
M_.param_names_long = char(M_.param_names_long, 'rhoepsKB');
M_.param_names = char(M_.param_names, 'BYSS_data');
M_.param_names_tex = char(M_.param_names_tex, 'BYSS\_data');
M_.param_names_long = char(M_.param_names_long, 'BYSS_data');
M_.param_names = char(M_.param_names, 'spread_data');
M_.param_names_tex = char(M_.param_names_tex, 'spread\_data');
M_.param_names_long = char(M_.param_names_long, 'spread_data');
M_.exo_det_nbr = 0;
M_.exo_nbr = 25;
M_.endo_nbr = 80;
M_.param_nbr = 83;
M_.orig_endo_nbr = 80;
M_.aux_vars = [];
M_.Sigma_e = zeros(25, 25);
M_.Correlation_matrix = eye(25, 25);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('FI_AINO16_rep_static');
erase_compiled_function('FI_AINO16_rep_dynamic');
M_.lead_lag_incidence = [
 0 48 0;
 1 49 128;
 2 50 129;
 3 51 130;
 4 52 0;
 5 53 0;
 6 54 0;
 7 55 131;
 0 56 0;
 0 57 0;
 0 58 0;
 8 59 132;
 9 60 133;
 0 61 0;
 10 62 0;
 11 63 134;
 0 64 0;
 0 65 0;
 12 66 0;
 0 67 0;
 13 68 0;
 14 69 0;
 0 70 135;
 0 71 0;
 0 72 0;
 0 73 0;
 0 74 0;
 15 75 136;
 16 76 0;
 17 77 137;
 0 78 0;
 0 79 0;
 0 80 0;
 18 81 0;
 0 82 0;
 19 83 0;
 0 84 0;
 0 85 138;
 0 86 0;
 20 87 0;
 0 88 0;
 0 89 0;
 0 90 0;
 21 91 139;
 22 92 0;
 23 93 0;
 24 94 0;
 25 95 0;
 26 96 0;
 27 97 0;
 28 98 0;
 29 99 0;
 30 100 140;
 31 101 0;
 32 102 0;
 33 103 0;
 34 104 0;
 35 105 0;
 0 106 0;
 36 107 0;
 37 108 0;
 38 109 0;
 0 110 0;
 39 111 0;
 40 112 0;
 0 113 0;
 41 114 141;
 42 115 0;
 43 116 0;
 44 117 0;
 0 118 0;
 0 119 0;
 0 120 0;
 45 121 0;
 0 122 0;
 46 123 0;
 0 124 0;
 47 125 0;
 0 126 0;
 0 127 0;]';
M_.nstatic = 31;
M_.nfwrd   = 2;
M_.npred   = 35;
M_.nboth   = 12;
M_.nsfwrd   = 14;
M_.nspred   = 47;
M_.ndynamic   = 49;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:25];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(80, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(25, 1);
M_.params = NaN(83, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 382;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
pp=load('all_parameters_February2016');
for i=1:M_.param_nbr;
name = deblank(M_.param_names(i,:));
if isfield(pp,name);
M_.params(i) = eval(['pp.' name]);
end;
end;
seepsZETAEUR=pp.seepsZETAEUR; 
seepsZETACH=pp.seepsZETACH;
seepsLAMW=pp.seepsLAMW;
seepsMU=pp.seepsMU;
seepsUPSILON=pp.seepsUPSILON;
seepsLAMBDALT=pp.seepsLAMBDALT;
seepsLAMBDAK=pp.seepsLAMBDAK;
seepsLAMBDACY=pp.seepsLAMBDACY;
seepsLAMBDACM=pp.seepsLAMBDACM;
seepsLAMBDAIY=pp.seepsLAMBDAIY;
seepsUPSILONMC=pp.seepsUPSILONMC;
seepsUPSILONX=pp.seepsUPSILONX;
seepsXX=pp.seepsXX;
seepsGF=pp.seepsGF;
seepsIG=pp.seepsIG;
seepshG=pp.seepshG;
seepsdS=pp.seepsdS;
seepsPIEW=pp.seepsPIEW;
seepsMW=pp.seepsMW;
seepsPOILS=pp.seepsPOILS;
seepsPRAWS=pp.seepsPRAWS;
seepsBankCapital=pp.seepsBankCapital;
seepsEPSB=pp.seepsEPSB;
seepsrEUR=pp.seepsrEUR; 
options_.solve_algo = 1;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
resid(1);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(2, 2) = (seepsLAMBDAK)^2;
M_.Sigma_e(25, 25) = (seepsrEUR)^2;
options_.irf = 20;
options_.nograph = 1;
options_.noprint = 1;
options_.order = 1;
var_list_=[];
var_list_ = 'y';
var_list_ = char(var_list_, 'cH');
var_list_ = char(var_list_, 'iH');
var_list_ = char(var_list_, 'x');
var_list_ = char(var_list_, 'm');
var_list_ = char(var_list_, 'tbY');
var_list_ = char(var_list_, 'rK');
var_list_ = char(var_list_, 'rwage');
var_list_ = char(var_list_, 'mcY');
var_list_ = char(var_list_, 'hF');
var_list_ = char(var_list_, 'ToT');
var_list_ = char(var_list_, 'pieY');
var_list_ = char(var_list_, 'rb');
var_list_ = char(var_list_, 'btot');
var_list_ = char(var_list_, 'bka');
var_list_ = char(var_list_, 'nwe');
var_list_ = char(var_list_, 'q');
var_list_ = char(var_list_, 'rs');
var_list_ = char(var_list_, 'rEUR');
info = stoch_simul(var_list_);
save('FI_AINO16_rep_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('FI_AINO16_rep_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('FI_AINO16_rep_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('FI_AINO16_rep_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('FI_AINO16_rep_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
