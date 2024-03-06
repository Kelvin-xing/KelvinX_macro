%%               Amsterdam Macroeconomics Summer School 2010
%                       Part II: Heterogeneous Agents
%                          University of Amsterdam
% 
%                            Thursday Assignment
%             Badly behaved higher-order pertubation solutions &
%                         why pruning is a bad idea
%--------------------------------------------------------------------------

%Cleaning
clear all; close all; clc

%% Parameters
%--------------------------------------------------------------------------

%Model parameters
beta  = 0.9;
gamma = 3;
mu    = 1.5;
sigma = 0.15;
r     = 0.03;

%Penalty parameters
penalty = '30';     %Type '10' for eta0 = 10 or '30' for eta0 = 30.
switch penalty
    case '10'
        eta0 = 10;  eta1 = 0.054;   eta2 = -0.0116;
    case '30'
        eta0 = 30;  eta1 = 0.045;   eta2 = -0.0018;
end

%Creating parameters.mod using system commands
delete parameters.mod                                                   %First delete old file, otherwise parameters are appended to old file.
system(sprintf('echo beta  = %0.14f; >> parameters.mod',beta )');       %"sprintf" creates text string. "system" executes text string as system command
system(sprintf('echo gamma = %0.14f; >> parameters.mod',gamma)');       %(as if typed in command window of operating system).
system(sprintf('echo mu    = %0.14f; >> parameters.mod',mu   )');
system(sprintf('echo sigma = %0.14f; >> parameters.mod',sigma)');
system(sprintf('echo r     = %0.14f; >> parameters.mod',r    )');
system(sprintf('echo eta0  = %0.14f; >> parameters.mod',eta0 )');
system(sprintf('echo eta1  = %0.14f; >> parameters.mod',eta1 )');
system(sprintf('echo eta2  = %0.14f; >> parameters.mod',eta2 )');

%% Steady state
%--------------------------------------------------------------------------

%Solving for steady state
init    = 0.05;
options = optimset('Display','Iter','MaxIter',1E5,'TolFun',1E-10,'TolX',1E-10);
yss     = mu;
ass     = fsolve(@(a) (a*r/(1+r)+yss)^-gamma*(1-(1+r)*beta)-(1+r)*(eta1*exp(-eta0*a)+eta2),init,options);
css     = ass*r/(1+r)+yss;

%Creating steadystate.mod using system commands
delete steadystate.mod                                             %First delete old file, otherwise parameters are appended to old file.
system(sprintf('echo y = %0.14f; >> steadystate.mod',yss)');       %"sprintf" creates text string. "system" executes text string as system command
system(sprintf('echo a = %0.14f; >> steadystate.mod',ass)');       %(as if typed in command window of operating system).
system(sprintf('echo c = %0.14f; >> steadystate.mod',css)');

%% Run Dynare
%--------------------------------------------------------------------------

dynare model.mod noclearall                 %Running Dynare without clearing memory
load dynarerocks                            %Loading decision rule using Wouter's alternative disp_dr.m,
                                            %so make sure you have replaced Dynare's standard disp_dr.m.

%Cleaning
!del model.m
!del model_static.m
!del model_dynamic.m
!del model_results.mat
!del model.log
!rmdir model

%% Policy function
%--------------------------------------------------------------------------

%Grid for cash on hand
xgrid = (1:0.01:2.5);
xhat = xgrid-ass-yss;

%QUESTION 1. PLOT THE POLICY FUNCTION. COMPLETE "ACHOICE = ..." AND YOU CAN
%PLOT FIGURE 1 BY UNCOMMENTING. NOTE THAT DYNARE RETURNS THE POLICY RULE FOR
%"A" AS FUNCTION OF "A(-1)" AND "E". HOW CAN YOU REWRITE THE POLICY RULE FOR
%"A" AS FUNCTION OF X INSTEAD? WHY IS THIS POSSIBLE?
%achoice = ...;

%Plot assets
% figure(1), plot(xgrid,achoice)
% title('Assets'), box off
% xlabel('\it{x_{t}}'), ylabel('\it{a_{t}}','Rotation',0)

%QUESTION 1. IN WHICH REGION IS E(x(t+1)) > x(t)? HOW CAN YOU VISUALIZE?
%hold all
%plot(...)

%% Simulation without pruning
%--------------------------------------------------------------------------

%Settings
simsize = 250;                                          %Number of periods

%Stochastics
randn('state',1038)                                     %State of generator
shock = sigma*randn(simsize,1);                         %Shock to dividend

%Simulation without pruning
ysim = mu+shock;
xsim = zeros(simsize,1);
asim = zeros(simsize,1);
for i=2:simsize,
    %QUESTION 2.1 AND 2.2. COMPLETE SIMULATION AND YOU CAN PLOT FIGURE 2
    %BY UNCOMMENTING.
    %xsim(i,1) = ...;
    %asim(i,1) = ...;
end

%Plot simulation
% figure(2)
% plot(1:simsize,asim)
% title('Assets'), box off
% xlabel('\it{t}'), ylabel('\it{a_{t}}','Rotation',0)

%% Simulation with pruning
%--------------------------------------------------------------------------

%Settings
simsize = 250;                                          %Number of periods

%Stochastics
randn('state',1038)                                     %State of generator
shock = sigma*randn(simsize,1);                         %Shock to dividend

%Simulation with pruning
ysim  = mu+shock;
xsim1 = zeros(simsize,1);   xsim2 = zeros(simsize,1);
asim1 = zeros(simsize,1);   asim2 = zeros(simsize,1);
for i = 2:simsize,
    %QUESTION 2.3. COMPLETE PRUNING SIMULATION AND YOU CAN PLOT FIGURE 3
    %BY UNCOMMENTING.
    %xsim1(i,1) = ...;
    %asim1(i,1) = ...;
    %xsim2(i,1) = ...;
    %asim2(i,1) = ...;    
end

%Plot simulation with pruning
% figure(3)
% plot(1:simsize,asim2)
% title('Assets pruned'), box off
% xlabel('\it{t}'), ylabel('\it{a_{t}}','Rotation',0)

%% Pruning policy "function"
%--------------------------------------------------------------------------

%QUESTION 2.4 AND 2.5. CREATE SCATTER PLOT WITH PRUNING POLICY "FUNCTION".
%ALSO INCLUDE REGULAR POLICY FUNCTION. HINT: TO GET PRUNING POLICY "FUNCTION"
%YOU SHOULD INCREASE THE NUMBER OF PERIODS IN YOUR SIMULATION TO SAY 25000. WHY?

%Scatter plot with pruning policy "function"
% figure(4)
% ...
% title('Pruning policy "function"')
% xlabel('\it{x_{t}}'), ylabel('\it{a_{t}}','Rotation',0)

%Include policy function
% hold all
% ...