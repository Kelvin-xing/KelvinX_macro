T = 1000;  % periods we wait for weights to converge

K = 5001;

agrid = nodeunif(K, smin(1), smax(1));


edges = [agrid(1); (agrid(2:end-1) + agrid(1:end-2))/2; agrid(end)];

state = gridmake(agrid, (1:1:k)');

N = size(state, 1); 
n = 1/N*ones(N,1);

% Decision rules

x   = funeval(cx, fspace, state); x = min(max(x, smin(1)),smax(1)); % don't allow outside bounds 

% Third, productive

xprime     = x;                      % add entrants
[~, xpbin] = histc(xprime, edges);   % bin by asset
xpbin      = [xpbin; (1:1:K)'];      % include all

% Iterate over measures

for t = 1 : T
    
nold = n;

nprime = nold;
Pprime = P(state(:,2),:);
Pprime = bsxfun(@times, nprime, Pprime); 

Pprime = [Pprime; zeros(K,k)];    % make sure all guys included

Pnew = zeros(K, k); 

for i = 1:k
    
   Pnew(:,i) = accumarray(xpbin,Pprime(:,i));    % add all guys in a given bin
    
end

n = reshape(Pnew, K*k, 1); 

%fprintf('%4i %6.2e \n',[t, norm(n-nold)]);

if norm(n-nold) < 1e-10, break, end

end

% Productive Sector

rr   = funeval(cr,fspace,state);
L   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(egrid(state(:,2)));
K   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(egrid(state(:,2)));
Y   = exp(egrid(state(:,2))).^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D   = Y - W.*L - (r+delta).*K;  
C   = D + (1+r)*state(:,1) - x; 

Ya = n'*Y;  % aggregate per capita output productive
Ka = n'*K; 
La = n'*L;

Aa = n'*state(:,1);  % assets in productive sector

TFP = Ya./(La.^alpha.*Ka.^(1-alpha)).^(eta);

TFPbest = (n'*exp(egrid(state(:,2)))).^(1-eta);




