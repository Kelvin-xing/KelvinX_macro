clear all;
clc;

% Adjust path to folder where replication file is stored
cd([cd '/EACZ_GEM03_rep']);

% Run replication dynare files

dynare EACZ_GEM03_rep;

ZZ_PIE4H = sqrt(oo_.var(1,1));
ZZ_DRNOMH = sqrt(oo_.var(2,2));
GDPGAPH = sqrt(oo_.var(3,3));
save EACZ_GEM03_rep_results.mat ZZ_PIE4H ZZ_DRNOMH GDPGAPH;
disp(' ');
disp(' ');
disp('Replication results:');
disp('Standard deviations of...');

disp(' ');
disp('Annual inflation');
disp('ZZ_PIE4H');
disp(ZZ_PIE4H);
disp(' ');
disp('First difference of interest rate');
disp('ZZ_DRNOMH');
disp(ZZ_DRNOMH);
disp(' ');
disp('Output gap');
disp('GDPGAPH');
disp(GDPGAPH);
