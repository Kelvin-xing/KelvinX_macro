% x=hermpol(ord,p);
% Purpose
% Create an ord-th order Hermit Polynomial 
% if p	is a vector then the polynomial will be created based 
%	on that vector
%
% ----------------------------------------------------------------------
function x=hermpol(ord,p)

    r	= size(p,1);
 	x	= ones(r,ord+1);	
	x(:,2)	= p;
	if ord >= 2;
        for i	= 3:ord+1;
           x(:,i)	= p.*x(:,i-1)-(i-1)*x(:,i-2);
        end;
	end;

% **********************************************************************
% **********************************************************************
