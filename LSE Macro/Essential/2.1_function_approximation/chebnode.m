% x=chebnode(n);
% Purpose
% Create n Chebyshev nodes
%
% November 9 1998
%
% ----------------------------------------------------------------------
function x=chebnode(n)

r	= max(1,n);
n	= (1:n)';
x	= cos( (pi*(n-0.5))/r );

% **********************************************************************

% **********************************************************************