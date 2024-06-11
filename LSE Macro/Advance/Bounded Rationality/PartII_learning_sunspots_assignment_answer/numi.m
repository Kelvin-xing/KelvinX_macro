% function yvec=numi(fun,par);
% Calculate the expected value of E[fun(epsi)] using Hermite Gaussian quadrature,
% when epsi is a normally distributed random variable with standard deviation sigma.
%
%

function yvec=numi(fun,par)


% Define the Nodes and weights for integration:

[nodes weight]  = hernodes(par.N_herm);

% -------------------------------------------------------

% Note that Hermite Gaussian quadrature uses the weighting function
% exp{-x^2}. Here, however, we have the density 
% (1./sigma*sqrt(2*pi))*exp{-y^2/2*sigma^2} and want to calculate the
% integral of fun(y)*(1./sigma*sqrt(2*pi))*exp{-y^2/2*sigma^2}
%
% Adjusting the quadrature procedure requires the following two steps.
% First a change of variables is introduced. In particular, y=sqrt(2)*sigma*x.
% This results in the weighting function (1/sigma*sqrt(2*pi))*exp{-x^2}. 
% Second, we have to multiply the integrand with the Jacobian of the 
% transformation. This gives the weighting function (1/sqrt(pi))*exp{-x^2}
% and calculating the integral of fun(sqrt(2)*sigma*x)*(1/sqrt(pi))*exp{-x^2}.
%

nodet   = par.vola*sqrt(2)*nodes;

temp = feval(fun,0,par);
vals = zeros(size(temp,1),par.N_herm);

% Evaluate the function at the transformed nodes

for hi  = 1:par.N_herm
	vals(:,hi)        = feval(fun,nodet(hi),par);
end

% Compute the integral
yvec    = 1/sqrt(pi) * vals*weight;

% **********************************************************************

% **********************************************************************
