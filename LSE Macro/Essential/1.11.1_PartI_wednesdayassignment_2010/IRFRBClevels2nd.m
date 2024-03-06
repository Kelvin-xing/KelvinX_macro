%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amsterdam Macroeconomics Summer School 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program computes IRFs for the second-order approximation. The
% purpose is to show that both initial values and future shock paths matter
% for the impulse responses when you are using the second order
% approximation
%

clc
clear

% These are the two parameters that will be loaded to Dynare files
nu    = 3;      
sigma = 0.007;  
save nuparam nu;     %This saves the value of nu to a file called nuparam
save sigparam sigma; %This saves the value of sigma to a file called sigparam
LengthIRF = 40;     %Sets the length of the IRF
T = LengthIRF+1;

% Generates a shock series. You will be asked to change the seed (666) to 
% something else and compare the impulse responses.  
randn('state',666);
shocks = sigma*randn(T,1);

shocks1 = [shocks(1,1)+sigma; shocks(2:end)];

% Solve the model with Dynare and load the decision matrix:
dynare RBClevels2nd.mod noclearall
load dynarerocks;

%--------------------------------------------------------------------------
% COMPUTING IRFs
%--------------------------------------------------------------------------

% Reserving space for variables that will hold the solution
k0   = zeros(T,1);
c0   = zeros(T,1);
z0   = zeros(T,1);

k1   = zeros(T,1);
c1   = zeros(T,1);
z1   = zeros(T,1);

% For the second-order approximation, the first row of the matrix "decision" 
% does not contain deterministic steady states anymore. Dynare prints the 
% correction in the second row of the decisio matrix. 
% This is why we have to correct for this:
c_ss = decision(1,1)-decision(2,1);
k_ss = decision(1,2)-decision(2,2);
z_ss = decision(1,3)-decision(2,3);

% To generate the impulse response, you could start at the steady state, 
% but also somewhere else. For the first order approximation, this will not
% matter, but for the second-order approximation, starting values are important.
% You should play around by changing these starting values, e.g., by
% starting 50% above the deterministic steady state for capital.
c0(1,1) = c_ss;
k0(1,1) = k_ss;
z0(1,1) = z_ss;

%This just ensures that both sets start at the same initial value
c1(1,1) = c0(1,1);
k1(1,1) = k0(1,1);
z1(1,1) = z0(1,1);

%The following recursion computes impulse responses for capital,
%consumption, and productivity. Note that the first element in the series
%is the steady state, which is why the recursion starts in period 2.

for i = 2:T
    %Vector S contains the "polynomial" that will be multiplied by the
    %coefficients given by Dynare's decision rules. The elements of "S" are
    %items listed on the left of policy and transition functions that
    %Dynare prints on screen.
    
    %This is the first set of series, which is computed only based on a
    %series of shocks. You can think of this as the "base case".
    S = [1,...
         0,...
         (k0(i-1)-k_ss(1,1)),...
         (z0(i-1)-z_ss(1,1)),...
         shocks(i-1),...
         (k0(i-1)-k_ss(1,1))^2,...
         (z0(i-1)-z_ss(1,1))*(k0(i-1)-k_ss(1,1)),...
         (z0(i-1)-z_ss(1,1))^2,...
         shocks(i-1)^2,...
         (k0(i-1)-k_ss(1,1))*shocks(i-1),...
         (z0(i-1)-z_ss(1,1))*shocks(i-1)];
    
    c0(i,1) = S*decision(:,1);
    k0(i,1) = S*decision(:,2);
    z0(i,1) = S*decision(:,3);
    
    %This is the second set of series, which is computed based on the same
    %shock series as the previous set, with the difference that a one
    %standard deviation shock has been added to the first element of the
    %shock series. You can think of this as the "shocked case".
    S = [1,...
         0,...
         (k1(i-1)-k_ss(1,1)),...
         (z1(i-1)-z_ss(1,1)),...
         shocks1(i-1),...
         (k1(i-1)-k_ss(1,1))^2,...
         (z1(i-1)-z_ss(1,1))*(k1(i-1)-k_ss(1,1)),...
         (z1(i-1)-z_ss(1,1))^2,...
         shocks1(i-1)^2,...
         (k1(i-1)-k_ss(1,1))*shocks1(i-1),...
         (z1(i-1)-z_ss(1,1))*shocks1(i-1)];
    
    c1(i,1) = S*decision(:,1);
    k1(i,1) = S*decision(:,2);
    z1(i,1) = S*decision(:,3);
end	

% The impulse-response is computed as the difference between the "base case" 
% and the "shocked case" 

c_irf = c1-c0;
k_irf = k1-k0;
z_irf = z1-z0;

figure
subplot(3,1,1),plot(c_irf(2:end,1)), xlim([1 40]), title('c')
subplot(3,1,2),plot(k_irf(2:end,1)), xlim([1 40]), title('k')
subplot(3,1,3),plot(z_irf(2:end,1)), xlim([1 40]), title('z')




