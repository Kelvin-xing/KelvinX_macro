%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************


clc
clear all

update = 0.5;
% the initial values used here are pretty bad so the value of the dampening
% coefficient cannot be too high (unless wcoef is low)

A_rational = 0; % = 1 rational o.w not

% true STRUCTURAL parameter values
alpha      = 0.7;
eta        = 0.1;
rho_z      = 0.8;
rho_i      = 0.2;
sig_z      = 0.01;
sig_i      = 0.1;
betta      = 0.99;
wcoef      = 0.5; % this corresponds to omega_1 in the program

save parametervalues alpha rho_z rho_i sig_z sig_i eta wcoef s_prob betta

%
% first run the representative agent version
%

dynare  heterobeliefs1repagent  noclearall
load    dynarerocks
disp('hit return to continue')
pause


% set initial values ALM as perceived by rational agent

coef_N_z1     = decision(6,3);
coef_N_z2     = decision(7,3);
coef_N_Nlag   = decision(5,3);

% using the following as initial values also works

% coef_N_z1     = 0;
% coef_N_z2     = 0;
% coef_N_Nlag   = 0;

% values perceived by type A firms

if A_rational == 1
    rho_z_A    = rho_z;
    rho_i_A    = rho_i;
    wcoef_A    = wcoef; 
    coef_N_z1_A   = coef_N_z1;
    coef_N_z2_A   = coef_N_z2;
    coef_N_Nlag_A = coef_N_Nlag;
    
else
    rho_z_A    = 0;
    rho_i_A    = 0;
    wcoef_A    = 0; 
    coef_N_z1_A   = 0.0;
    coef_N_z2_A   = 0.0;
    coef_N_Nlag_A = 0.0;
    
end

save parametervalues alpha rho_z rho_z_A rho_i rho_i_A sig_z sig_i eta wcoef wcoef_A s_prob betta

% coefficients of wage expectation

T  = 10000; %number of observations in time series
T1 = 1001;  %number of observations discarded
I = 1000;   %number of agents



% set up memory
N_simul  = zeros(T,1);
z_values = zeros(T,2);

% generate productivity values

randn('seed',20100814)
shocks = randn(T,2);
seed_save = randn('seed'); % here we check what the seed is at this point so
                           % so we always continue the one draw we started 
    

z_values(1,:) = zeros(1,2);
for t = 2:T;
    z_values(t,:) = rho_z*z_values(t-1,:) + sig_z*shocks(t,:);
    %!!! these are the values relative to steady state
end

% set up part of the explanatory variables that is fixed

X_fixed = z_values;

%
% start the main iteration loop
%

error = 100;

while error > 0.0001;
    
randn('seed',seed_save)  % in every loop we use the same seed and thus the 
                         % same random numbers
    
%
% get policy functions with Dynare
%

save    aggregatelaw_N  coef_N_z1 coef_N_z2 coef_N_Nlag coef_N_z1_A coef_N_z2_A coef_N_Nlag_A 
                        
dynare  heterobeliefs1  noclearall
load    dynarerocks


% initialize the cross section in period t=1
% we assume that initially everybody is the same which is done for
% convenience (but means we clearly have to discard some initial values)

N_simul(1) = 1;
n_1_R_old = 0.5*ones(I,1);
n_2_R_old = 0.5*ones(I,1);
n_1_A_old = 0.5*ones(I,1);
n_2_A_old = 0.5*ones(I,1);
z_i_R_old =    zeros(I,1);
z_i_A_old =    zeros(I,1);


for t = 2:T
    shocks_R = randn(I,1);
    shocks_A = randn(I,1);
    z_i_R    = rho_i*z_i_R_old + sig_i*shocks_R;
    z_i_A    = rho_i*z_i_A_old + sig_i*shocks_A;
    %!!! these are relative to steady state
    
    n_new =          [ z_i_R*decision(12,1)  z_i_R*decision(12,2)  z_i_A*decision(12,3)  z_i_A*decision(12,4)];
    n_new = n_new +  [n_1_R_old-0.5 n_2_R_old-0.5 n_1_A_old-0.5 n_2_A_old-0.5]*decision(6:9,1:4) ;
    n_new = n_new +  ...
                  +  [ones(I,1)*decision(1,1) ones(I,1)*decision(1,2) ones(I,1)*decision(1,3) ones(I,1)*decision(1,4)];
    n_new = n_new +  (N_simul(t-1)-1)* ...
                     [ones(I,1)*decision(2,1) ones(I,1)*decision(2,2) ones(I,1)*decision(2,3) ones(I,1)*decision(2,4)];
    n_new = n_new +  z_values(t,1)* ...
                     [ones(I,1)*decision(10,1) ones(I,1)*decision(10,2) ones(I,1)*decision(10,3) ones(I,1)*decision(10,4)];
    n_new = n_new +   z_values(t,2)* ...
                     [ones(I,1)*decision(11,1) ones(I,1)*decision(11,2) ones(I,1)*decision(11,3) ones(I,1)*decision(11,4)];
                              
    N_simul(t) = sum(sum(n_new))/(I*2);
    %disp(N_simul(t));
    %pause
    z_i_R_old = z_i_R;
    z_i_A_old = z_i_A;
    n_1_R_old = n_new(:,1);
    n_2_R_old = n_new(:,2);
    n_1_A_old = n_new(:,3);
    n_2_A_old = n_new(:,4);
end    

N_simul = N_simul - mean(N_simul);
X = [X_fixed(T1:end,:) N_simul(T1-1:end-1)];
Y = [N_simul(T1:end)];
coef=inv(X'*X)*X'*Y;
disp(coef');
coef_old = [coef_N_z1 coef_N_z2 coef_N_Nlag]';
disp(coef_old')
coef_N_z1   = update*coef(1) + (1-update)*coef_N_z1;
coef_N_z2   = update*coef(2) + (1-update)*coef_N_z2;
coef_N_Nlag = update*coef(3) + (1-update)*coef_N_Nlag;


if A_rational == 1
    coef_N_z1_A   = coef_N_z1;
    coef_N_z2_A   = coef_N_z2;
    coef_N_Nlag_A = coef_N_Nlag;
else
    coef_N_z1_A   = 0.0;
    coef_N_z2_A   = 0.0;
    coef_N_Nlag_A = 0.0;
end

error = sum(abs(coef-coef_old));
disp(error)
pause(3)

end