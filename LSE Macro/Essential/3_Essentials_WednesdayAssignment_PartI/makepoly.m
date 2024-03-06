% xmat=makepoly(ord,x);
%
%            For the (Tx2) matrix x and the (1x2) vector ord
%            this creates a [T x ord(1)*ord(2)] matrix where 
%            each row contains the ord(1)*ord(2) bivariate Hermite 
%            polynomial terms.
%
% ------------------------------------------------------------------

function xmat=makepoly(ord,x)

po_k    = ord(1);
po_t    = ord(2);
xk      = hermpol(po_k,x(:,1));
xt      = hermpol(po_t,x(:,2));
xmat    = zeros(size(x,1),(po_k+1)*(po_t+1));
for kj  = 0:po_k;
	for tj  = 0:po_t;
		xmat(:,kj*(po_t+1)+tj+1) = xk(:,kj+1).*xt(:,tj+1);
	end;
end;

% **********************************************************************

% **********************************************************************
