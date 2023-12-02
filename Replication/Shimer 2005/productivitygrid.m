function grid = productivitygrid(mu, sigma)
% Size of the grid
n = round(mu/sigma) - 1; % one side grid

% Create the grid
p = lognrnd(mu,sigma,1,2*n+1);
grid = sort(p);