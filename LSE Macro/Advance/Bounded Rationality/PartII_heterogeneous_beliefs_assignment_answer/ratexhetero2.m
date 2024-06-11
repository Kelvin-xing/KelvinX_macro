%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
c
clc
clear all

update = 0.9;

A_rational = 1; % = 1 rational o.w not
% although the program still calculates the expectations of type A agents
% their expectations do not matter since their 
% policy rules are not dependent on expectations

% true STRUCTURAL parameter values
alpha      = 0.7;
eta        = 0.1;
rho_z      = 0.8;
rho_i      = 0.2;
sig_z      = 0.01;
sig_i      = 0.1;
betta      = 0.99;
s_prob     = 0.02;
wcoef      = 0.1; % this corresponds to omega_1 in the program

% values perceived by type A firms

if A_rational == 1
    rho_z_A    = 0.8;
    rho_i_A    = 0.2;
    wcoef_A    = 0.1; 
else
    rho_z_A    = 0;
    rho_i_A    = 0;
    wcoef_A    = 0; 
end

% coefficients of wage expectation

save parametervalues alpha rho_z rho_z_A rho_i rho_i_A sig_z sig_i eta wcoef wcoef_A s_prob betta
T  = 5000; %number of observations in time series
T1 = 501;  %number of observations discarded
I = 500;   %number of agents of each type
s_num = s_prob*2*I; %number of switchers must be even number !!!!

% set initial values ALM as perceived by rational agent
coef_N_z1     = 0.0001;
coef_N_z2     = 0.0001;
coef_N_Nlag   = 0.0001;
% values ALM as perceived by irrational agent
%    these could be initial values and updated
%    or they could be kept fixed.
coef_N_z1_A   = 0.0;
coef_N_z2_A   = 0.0;
coef_N_Nlag_A = 0.0;

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
dynare  heterobeliefs2  noclearall
load    dynarerocks

% initialize the cross section in period t=1
% we assume that initially everybody is the same which is done for
% convenience (but means we clearly have to discard some initial values)

N_simul(1) = 1;
n_1_old = 0.5*ones(2*I,1);
n_2_old = 0.5*ones(2*I,1);
z_i_old =    zeros(2*I,1);
n_new   =    zeros(2*I,2);

R_indic = [ ones(I,1); -ones(I,1)]; % first I agents are rational other not

for t = 2:T
    shocks = randn(2*I,1);
    z_i    = rho_i*z_i_old + sig_i*shocks;    
    %!!! these are relative to steady state
    
    for i = 1:2*I
        if R_indic(i) > 0
            n_new(i,:) =  [ z_i(i)*decision(10,1)  z_i(i)*decision(10,2) ];
            n_new(i,:) =  n_new(i,:) + [n_1_old(i)-0.5 n_2_old(i)-0.5 ]*decision(6:7,1:2) ;
            n_new(i,:) =  n_new(i,:) +  ...
                       +  [decision(1,1) decision(1,2) ];
            n_new(i,:) =  n_new(i,:) + (N_simul(t-1)-1)* ...
                          [decision(2,1) decision(2,2) ];
            n_new(i,:) =  n_new(i,:) +  z_values(t,1)* ...
                          [decision(8,1) decision(8,2) ];
            n_new(i,:) =  n_new(i,:) +   z_values(t,2)* ...
                          [decision(9,1) decision(9,2) ];
        else
            n_new(i,:) =  [decision(1,1) decision(1,2) ];
        end

    end
    
    z_i_old = z_i;
    n_1_old = n_new(:,1);
    n_2_old = n_new(:,2);    
    
    N_simul(t) = sum(sum(n_new))/(I*2);
    %disp(N_simul(t));
    %pause


    % determine which agents are going to switch
    count_pos = 0;
    count_neg = 0;
    count     = 1;
    rand('seed',20100814)
%   an earlier version of the program did not set the seed. But 
%   randperm() calls rand() and you have to set separate seeds for 
%   rand and randn   
    
    ss  = randperm(2*I)';
    while count_pos+count_neg < s_num
        if R_indic(ss(count))>0 
            if count_pos < s_num/2 
                R_indic(ss(count))=-R_indic(ss(count));
                count_pos = count_pos+1;                
            end
        else
            if count_neg < s_num/2
                R_indic(ss(count))=-R_indic(ss(count));
                count_neg = count_neg+1;
            end
        end

        count = count+1;

    end
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
%pause

end