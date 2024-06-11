% function yvec=numi(fun,par);
% Calculate the expected value of E[fun(epsi)] 
% when epsi is has discrete support (+1 & -1)
%
%

function yvec=numi_case1(fun,par)


% Define the Nodes and weights for integration:

nodet = [-1;1];
weight = [0.5;0.5];

temp = feval(fun,0,par);
vals = zeros(size(temp,1),2);

% Evaluate the function at the transformed nodes

for hi  = 1:2
	vals(:,hi)        = feval(fun,nodet(hi),par);
end

% Compute the integral
yvec    = vals*weight;

% **********************************************************************

% **********************************************************************
