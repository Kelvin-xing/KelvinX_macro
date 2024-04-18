T = 1000;  % periods we wait for weights to converge
N = 5001;

curv   = 0.25;
agrid  = nodeunif(N, smin(1).^curv, smax(1).^curv).^(1/curv);

agrid(1)   = smin(1);  % may not be exactly equal because of numerical issues
agrid(end) = smax(1); 

edges  = [agrid(1); (agrid(2:end-1) + agrid(1:end-2))/2; agrid(end)];

state  = gridmake(agrid,  znode);
stated = gridmake(agrid,  (1:1:k)');

Ns  = size(state, 1); 

nm = 0.70*1/Ns*ones(Ns,1);
nt = 0.30*1/Ns*ones(Ns,1);

% Decision rules

xm   = funeval(cx(:,1), fspace, state); xm = min(max(xm, smin(1)), smax(1)); % don't allow outside bounds 
xt   = funeval(cx(:,2), fspace, state); xt = min(max(xt, smin(1)), smax(1)); % don't allow outside bounds 


Phi  = funbas(fspace,  [xt, state(:,2)]);
st   = Phi*cst(:,1) > Phi*cst(:,2) & xt  - kappa*exp(state(:,2)) > smin(1);   % switching decision of traditional sector producer     
xt(st==1) = xt(st==1) - kappa*exp(state(st==1,2));             

% Sort savings into bins

% First, traditional: a) trad. who stay, b) mod who join, c) new guys

xtprime     = [xt(st==0); xm; smin(1)];    
[~, xtpbin] = histc(xtprime, edges);    % bin by asset
xtpbin      = [xtpbin; (1:1:N)'];       % include all N

% Second, modern

xmprime     = [xm; xt(st==1)];          % add entrants
[~, xmpbin] = histc(xmprime, edges);    % bin by asset
xmpbin  = [xmpbin; (1:1:N)'];           % include all

% Iterate over measures

P1 = eye(k);                % transition probability for those who keep productivity
P2 = repmat(wnode',k,1);    % transition probability for those who switch productivity: do it like this if later want some memory
P3 = rhot*P1 + (1-rhot)*P2; % transition probability of traditional who switch 
P4 = rhom*P1 + (1-rhom)*P2; % transition probability of traditional who stay

for t = 1 : T
    
    % first, traditional
    
ntold = nt;
nmold = nm;

ntprime     = [1/mu*ntold(st==0);       1/mu*(1-rhom)*nmold;  1-1/mu];
Ptprime     = [P4(stated(st==0,2),:);   P2(stated(:,2),:);    wnode'];
Ptprime     = bsxfun(@times, ntprime,   Ptprime); 

Ptprime = [Ptprime; zeros(N,k)];    % make sure all guys included

Ptnew = zeros(N, k); 

for i = 1:k
    
   Ptnew(:,i) = accumarray(xtpbin, Ptprime(:,i));    % add all guys in a given bin
    
end

nt = reshape(Ptnew, N*k, 1); 

    % second, modern

nmprime = [1/mu*rhom*nmold;     1/mu*ntold(st==1)];
Pmprime = [P1(stated(:,2),:);   P3(stated(st==1,2),:)];
Pmprime = bsxfun(@times, nmprime, Pmprime); 

Pmprime = [Pmprime; zeros(N,k)];    % make sure all guys included

Pmnew = zeros(N, k); 

for i = 1:k
    
   Pmnew(:,i) = accumarray(xmpbin,Pmprime(:,i));    % add all guys in a given bin
    
end

nm = reshape(Pmnew, N*k, 1); 

%fprintf('%4i %6.2e %6.2e  \n',[t, norm(nt-ntold), norm(nm-nmold)]);

if norm(nt-ntold) < 1e-9 && norm(nm-nmold) < 1e-9, break, end

end

% Modern Sector

rr   = funeval(cr, fspace, state);
Lm   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(state(:,2) + phi);
Km   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(state(:,2) + phi);
Ym   = exp(state(:,2) + phi).^(1-eta).*(Lm.^alpha.*Km.^(1-alpha)).^eta;

Dm   = Ym - W.*Lm - (r+delta).*Km;  

Cm   = Dm + (1+r)*state(:,1) - xm; 

% Traditional Sector

Lt  = (eta/W)^(1/(1-eta))*exp(state(:,2));
Yt  = (eta/W)^(eta/(1-eta))*exp(state(:,2));
Dt  = Yt - W*Lt; 
Ct  = Dt + (1+r)*state(:,1) - (xt + kappa.*exp(state(:,2)).*(st==1));   % we subtracted kappa from these guys, so add it back

Yma = nm'*Ym;  % aggregate per capita output modern
Kma = nm'*Km;  
Lma = nm'*Lm;

Yta = nt'*Yt;  % traditional
Lta = nt'*Lt;

Ama = nm'*state(:,1);  % assets in modern sector
Ata = nt'*state(:,1);  % assets in traditional sector

TFPm = Yma./(Lma.^alpha.*Kma.^(1-alpha)).^(eta);
TFPt = Yta./Lta.^eta;

TFPmbest = (nm'*exp(state(:,2) + phi)).^(1-eta);

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Kma,        Ama + Ata]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [Lma + Lta,  Lbar]);
fprintf('\n');

D    = Kma - (Ama + Ata);    % Net Foreign Debt Position
Ya   = Yma + Yta; 
Ca   = Ya + (1 - delta - mu)*Kma - nt'*((st==1).*exp(state(:,2)))*kappa  - (1 + r - mu)*D; 
Ca2  = nt'*Ct + nm'*Cm + W*Lbar;   % add consumption of all agents

sdY = sqrt(1/sum(nm)*nm'*(log(Ym)-1/sum(nm)*nm'*log(Ym)).^2);

data = sortrows([nm/sum(nm), Ym], 2); 
cweight = cumsum(data(:,1)); 
f10 = data(cweight>.90,1)'*data(cweight>.9,2)/(data(:,1)'*data(:,2));

emplpremium = 1/sum(nm)*nm'*Lm/(1/sum(nt)*nt'*Lt);

fprintf('\n');
disp('Model Implications')
fprintf('\n');
fprintf('Fraction modern               = %9.3f \n',          sum(nm));
fprintf('Fraction output in modern     = %9.3f \n',          Yma/(Yma + Yta));
fprintf('Fraction labor in modern      = %9.3f \n',          Lma/(Lma + Lta));
fprintf('\n');

fprintf('Fraction new entrants         = %9.3f \n',          1/mu*sum(nt(:).*st(:))/sum(nm(:)));
fprintf('Std. dev. Y                   = %9.3f \n',          sdY);
fprintf('Fraction revenue top 10 perc  = %9.3f \n',          f10);
fprintf('Modern to Trad L ratio        = %9.2f  %9.2f \n',   [emplpremium, 5])

fprintf('\n');

finvest = nt'*((st==1).*exp(state(:,2)))*kappa;

fprintf('Output                        = %9.3f \n',          Ya/(Lta+Lma));   % per worker in case L doesnt clear
fprintf('Consumption                   = %9.3f \n',          Ca/(Lta+Lma));
fprintf('Investment                    = %9.3f \n',          ((delta+mu-1)*Kma + finvest)/(Lta+Lma));
fprintf('Fraction Fixed Investment     = %9.3f \n',          finvest/ ((delta+mu-1)*Kma + finvest));
fprintf('TFP Traditional               = %9.3f \n',          TFPt);
fprintf('TFP Modern                    = %9.3f \n',          TFPm);
fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPmbest/TFPm)*100);
