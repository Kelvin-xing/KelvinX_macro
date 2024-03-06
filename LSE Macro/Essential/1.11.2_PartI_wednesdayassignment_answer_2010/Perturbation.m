%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amsterdam Macroeconomics Summer School 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program solves the standard RBC model in four different ways, using
% perturbation methods. The four ways are the following:
%       - model in levels:
%               - 1st order
%               - 2nd order
%       - model in logs:
%               - 1st order
%               - 2nd order
% 
% For each solution it simulates capital and consumption series and plots them

clc
clear

% These are the two parameters that will be loaded to Dynare files
nu    = 3;      
sigma = 0.007;  
save nuparam nu;     %This saves the value of nu to a file called nuparam
save sigparam sigma; %This saves the value of sigma to a file called sigparam


% For each solution, the program will simulate artificial time series for
% capital and consumption and store them. 

T=1000;                 %Length of the simulated series
capstore=zeros(T,4);    %Storing series for capital
constore=zeros(T,4);    %Storing series for consumption

% The shock series will be the same for each solution. 
randn('state',666);
shocks = sigma*randn(T,1);

% Reserving space for variables that will hold the solution
k    = zeros(T,1);
c    = zeros(T,1);
z   = zeros(T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELS IN LEVELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST ORDER:

dynare RBClevels1st.mod noclearall
% The command noclearall means that variables defined earlier (like T) are not deleted 

% Now load the decision rules given by Dynare. This command loads the matrix 
% "decision" with coefficients of policy functions computed by Dynare
% (note: this requires the updated disp_dr.m file)
load dynarerocks

% The first row of the matrix "decision" contains steady states:
c_ss = decision(1,1);
k_ss = decision(1,2);
z_ss = decision(1,3);

% To generate the artificial data, we will start at the steady state.
c(1,1) = c_ss;
k(1,1) = k_ss;
z(1,1) = z_ss;

%The following recursion computes artificial time series for capital,
%consumption, and productivity. Note that the first element in the series
%is the steady state, which is why the recursion starts in period 2.

for i = 2:T
    %Vector S contains the "polynomial" that will be multiplied by the
    %coefficients given by Dynare's decision rules. The elements of "S" are
    %items listed on the left of policy and transition functions that
    %Dynare prints on screen.
    S = [1, k(i-1)-k_ss, z(i-1)-z_ss, shocks(i)];
    
    %Here the multiplication with coefficients is performed:
    c(i,1) = S*decision(:,1);
    k(i,1) = S*decision(:,2);
    z(i,1) = S*decision(:,3);
    
end	

% The following stores the artificial series generated using
% first-order approximation:
capstore(:,1)=k;
constore(:,1)=c;

%--------------------------------------------------------------------------
% SECOND ORDER

dynare RBClevels2nd.mod noclearall
% The command noclearall means that variables defined earlier (like T) are not deleted 

% Now load the decision rules given by Dynare. This command loads the matrix 
% "decision" with coefficients of policy functions computed by Dynare
% (note: this requires the updated disp_dr.m file)
load dynarerocks

% For the second-order approximation, the first row of the matrix "decision" 
% does not contain deterministic steady states anymore. Dynare prints the 
% correction in the second row of the decisio matrix. 
% This is why we have to correct for this:
c_ss = decision(1,1)-decision(2,1);
k_ss = decision(1,2)-decision(2,2);
z_ss = decision(1,3)-decision(2,3);

% To generate the artificial data, we will start at the deterministic steady state.
% Note that for the second-order approximation, starting values are important!
c(1,1) = c_ss;
k(1,1) = k_ss;
z(1,1) = z_ss;

%Below you should write the recursion that computes artificial time series 
%for capital, consumption, and productivity using Dynare's decision rules.
%Note that the first element in the series is the steady state, which is 
%why the recursion starts in period 2.

for i = 2:T
    
end	

% The following stores the artificial series generated using
% second-order approximation:
capstore(:,2)=k;
constore(:,2)=c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELS IN LOGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST ORDER

% Use Dynare to solve here the model you have rewritten in logs.
% (Don't forget to add noclearall command.)



% Now load the decision rules given by Dynare. This command loads the matrix 
% "decision" with coefficients of policy functions computed by Dynare
% (note: this requires the updated disp_dr.m file)
load dynarerocks

% The first row of the matrix "decision" contains steady states:
c_ss = decision(1,1);
k_ss = decision(1,2);
z_ss = decision(1,3);

% To generate the artificial data, we will start at the steady state.
c(1,1) = c_ss;
k(1,1) = k_ss;
z(1,1) = z_ss;

%The following recursion computes artificial time series for capital,
%consumption, and productivity. Note that the first element in the series
%is the steady state, which is why the recursion starts in period 2.

for i = 2:T
    %Vector S contains the "polynomial" that will be multiplied by the
    %coefficients given by Dynare's decision rules. The elements of "S" are
    %items listed on the left of policy and transition functions that
    %Dynare prints on screen.
    S = [1, k(i-1)-k_ss, z(i-1)-z_ss, shocks(i)];
    
    %Here the multiplication with coefficients is performed:
    c(i,1) = S*decision(:,1);
    k(i,1) = S*decision(:,2);
    z(i,1) = S*decision(:,3);
    
end	

% The following stores the artificial series generated using
% first-order approximation. Note that because the model is now in logs, we
% have to exp() the series to obtain the levels.
capstore(:,3)=exp(k);
constore(:,3)=exp(c);

%--------------------------------------------------------------------------
% SECOND ORDER

% Use Dynare to solve here the model you have rewritten in logs. 
% (Don't forget to add noclearall command.) 




% Now load the decision rules given by Dynare. This command loads the matrix 
% "decision" with coefficients of policy functions computed by Dynare
% (note: this requires the updated disp_dr.m file)
load dynarerocks

% For the second-order approximation, the first row of the matrix "decision" 
% does not contain deterministic steady states anymore. Dynare prints the 
% correction in the second row of the decisio matrix. 
% This is why we have to correct for this:
c_ss = decision(1,1)-decision(2,1);
k_ss = decision(1,2)-decision(2,2);
z_ss = decision(1,3)-decision(2,3);

% To generate the artificial data, we will start at the deterministic steady state.
% Note that for the second-order approximation, starting values are important!
c(1,1) = c_ss;
k(1,1) = k_ss;
z(1,1) = z_ss;

%Below you should write the recursion that computes artificial time series 
%for capital, consumption, and productivity using Dynare's decision rules.
%You could actually copy what you have written for the model in levels.
%Note that the first element in the series is the steady state, which is 
%why the recursion starts in period 2.

for i = 2:T

end	

% The following stores the artificial series generated using
% second-order approximation. Note that because the model is now in logs, we
% have to exp() the series to obtain the levels.
capstore(:,4)=exp(k);
constore(:,4)=exp(c);


% figure
% plot(capstore)
% legend('1st order levels','2nd order levels','1st order logs','2nd order logs',0)
% title('Artificial series for capital')
% figure
% plot(constore)
% legend('1st order levels','2nd order levels','1st order logs','2nd order logs',0)
% title('Artificial series for consumption')