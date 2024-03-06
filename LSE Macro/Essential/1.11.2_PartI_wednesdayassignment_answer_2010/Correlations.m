%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amsterdam Macroeconomics Summer School 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this program is to investigate the effects of different
% shocks on correlations of some key business cycle statistics

clc
clear

T=1000;                 %Length of the simulated series
sigma = 0.007;

% The shock series will be the same for each case. 
randn('state',666);
shocks = sigma*randn(T,1);

% Reserving space for variables that will hold the solution
k    = zeros(T,1);
c    = zeros(T,1);
inv  = zeros(T,1);
y    = zeros(T,1);
z    = zeros(T,1);

outstore = zeros(T,2);    %Storing series for output
constore = zeros(T,2);    %Storing series for consumption
invstore = zeros(T,2);    %storing series for investment


dynare RBCTechnology.mod noclearall

load dynarerocks

% The first row of the matrix "decision" contains steady states:
c_ss   = decision(1,1);
k_ss   = decision(1,2);
inv_ss = decision(1,3);
y_ss   = decision(1,4);
z_ss   = decision(1,5);

% Start at the steady state.
c(1,1)   = c_ss;
k(1,1)   = k_ss;
inv(1,1) = inv_ss;
y(1,1)   = y_ss;
z(1,1)   = z_ss;

% The loop that generates artificial data series 
for i = 2:T
    S = [1, k(i-1)-k_ss, z(i-1)-z_ss, shocks(i-1)];
    
    c(i,1)   = S*decision(:,1);
    k(i,1)   = S*decision(:,2);
    inv(i,1) = S*decision(:,3);
    y(i,1)   = S*decision(:,4);
    z(i,1)   = S*decision(:,5);
end	

%Store the computed series in the first column:
outstore(:,1) = y;
constore(:,1) = c;
invstore(:,1) = inv;

% Now repeat the exercise that was done above with the model with 
% investment technology shocks. You have to do the following: Solve
% RBCInvestment with Dynare (do not forget the noclearall command), 
% load the coefficient matrix (decision rules) and write a loop that
% generates artificial time series for y, c, and inv.

dynare ...








% This stores the artificial series you just computed in the second column:
outstore(:,2) = y;
constore(:,2) = c;
invstore(:,2) = inv;

% Business cycle statistics:

HPyTech   = outstore(:,1)-hpfilter(outstore(:,1),1600);
HPcTech   = constore(:,1)-hpfilter(constore(:,1),1600);
HPiTech   = invstore(:,1)-hpfilter(invstore(:,1),1600);

HPyInvest = outstore(:,2)-hpfilter(outstore(:,2),1600);
HPcInvest = constore(:,2)-hpfilter(constore(:,2),1600);
HPiInvest = invstore(:,2)-hpfilter(invstore(:,2),1600);

% Correlations

CorrTech = corr([HPyTech,HPcTech,HPiTech]);
disp('After a TFP shock we obtain the following correlations:')
disp('Output  Consumption  Investment')
disp(CorrTech)


CorrInvest = corr([HPyInvest,HPcInvest,HPiInvest]);
disp('After an investment technology shock we obtain the following correlations:')
disp('Output  Consumption  Investment')
disp(CorrInvest)