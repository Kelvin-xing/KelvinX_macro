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

L = eta^(1/(1-eta))*(W/alpha)^(-gamma)*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta))).*exp(egrid(state(:,2)));
K = eta^(1/(1-eta))*((rr+delta)/(1-alpha)).^(-gamma).*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma*(rr+delta).^(1-gamma)).^(1/(1-gamma)*(gamma-1/(1-eta))).*exp(egrid(state(:,2)));
Y = eta^(eta/(1-eta))*(alpha^gamma*W^(1-gamma) + (1-alpha)^gamma*(rr+delta).^(1-gamma)).^(1/(1-gamma)*eta/(eta-1)).*exp(egrid(state(:,2)));  

D   = Y - W.*L - (r+delta).*K;  
C   = D + (1+r)*state(:,1) - x; 

Ya = n'*Y;  % aggregate per capita output productive
Ka = n'*K; 
La = n'*L;

Aa = n'*state(:,1);  % assets in productive sector


disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Ka,  Aa + Aw]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [La,  Lbar]);
fprintf('\n');



D    = Ka - (Aa + Aw);                 % Net Foreign Debt Position
Ca   = Ya - delta*Ka - r*D; 
Ca2  = n'*C + W*Lbar + r*Aw;   % add consumption of all agents

fprintf('\n');
disp('Model Implications')

fprintf('Output                        = %9.3f \n',          Ya/La);   % per worker in case L doesnt clear
fprintf('Consumption                   = %9.3f \n',          Ca/La);
fprintf('Investment                    = %9.3f \n',          delta*Ka/La);
