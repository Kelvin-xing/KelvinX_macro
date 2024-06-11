%
% very simple program to solve the standard growth model with projection
% methods
%
% simple grid (although ideally Chebyshev nodes should be used)
% regular polynomial (although ideally Cheb orth. pol should be used)
%
% the program uses two functions
% consfun(k,z,coef) evaluates consumption
% griderror(coef,k_grid,z_grid) calculates sum of squared Euler eq errors 

clear
clc
global beta alpha depr gamma rho sigma q_nodes q_weights k_grid z_grid

% parameter values

beta  = 0.99;
alpha = 0.33;
depr  = 0.025;
gamma = 4;

rho   = 0.95;
sigma = 0.01;

% calculate steady state

k_ss = (beta*alpha/(1-beta*(1-depr)))^(1/(1-alpha));

% define grid

k_low    =  0.5*k_ss;
k_high   =  1.5*k_ss;
z_low    = -3*sqrt(sigma^2/(1-rho^2));
z_high   =  3*sqrt(sigma^2/(1-rho^2));

k_number =  10;
z_number =  10;

k_step   =  (k_high-k_low)/(k_number-1);
z_step   =  (z_high-z_low)/(z_number-1);

k_grid  = (k_low:k_step:k_high)';
z_grid  = (z_low:z_step:z_high)';
z_grid  = exp(z_grid);

% generate nodes and weights for numerical integration

q_number = 5;
[q_nodes,q_weights]= hernodes(q_number);

% initial value

coef_in = [log(1-alpha*beta); alpha; 1; 0; 0; 0];

% find solution by minimizing the errors on the grid

options = optimset('MaxFunEvals',100000,'MaxIter',1000000,'TolFun',0.0000001,'TolX',0.0000001);
coef_out = fminsearch(@griderror,coef_in,options);

% plot the consumption choice as a function of k (for 3 values of z)

k_grid_fine = (k_low:(k_high-k_low)/1000:k_high)';
plot(k_grid_fine,consfun(k_grid_fine,0.9,coef_out), ...
     k_grid_fine,consfun(k_grid_fine,1.0,coef_out), ...
     k_grid_fine,consfun(k_grid_fine,1.1,coef_out))