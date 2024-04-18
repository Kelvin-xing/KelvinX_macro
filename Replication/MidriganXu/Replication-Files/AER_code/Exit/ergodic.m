T = 1000;  % periods we wait for weights to converge

N = 5001;

agridp = nodeunif(N, sminp(1), smaxp(1));
agridu = nodeunif(N, sminu(1), smaxu(1));
agridt = nodeunif(N, smint(1), smaxt(1));

edgesp = [agridp(1); (agridp(2:end-1) + agridp(1:end-2))/2; agridp(end)];
edgesu = [agridu(1); (agridu(2:end-1) + agridu(1:end-2))/2; agridu(end)];
edgest = [agridt(1); (agridt(2:end-1) + agridt(1:end-2))/2; agridt(end)];

statep = gridmake(agridp, (1:1:k)');
stateu = gridmake(agridu, (1:1:k)');
statet = gridmake(agridt, (1:1:k)');

Np = size(statep, 1); 
Nu = size(stateu, 1); 
Nt = size(statet, 1);

np = 0.25*1/Np*ones(Np,1);
nu = 0.50*1/Nu*ones(Nu,1);
nt = 0.25*1/Nt*ones(Nt,1);

% Decision rules

xp   = funeval(cxp, fspacep, statep); xp = min(max(xp, sminp(1)),smaxp(1)); % don't allow outside bounds 
xu   = funeval(cxu, fspaceu, stateu); xu = min(max(xu, sminu(1)),smaxu(1)); % don't allow outside bounds 
xt   = funeval(cxt, fspacet, statet); xt = min(max(xt, smint(1)),smaxt(1)); % don't allow outside bounds 


Phi  = funbas(fspaceu, [xu, stateu(:,2)]);
su = Phi*csu(:,1) > Phi*csu(:,2) & xu - kappap > sminp(1);           % switching decision of unproductive producer       
xu(su==1) = xu(su==1) - kappap;


Phi  = funbas(fspacet,        [xt, statet(:,2)]);
inj  = funeval(cinj, fspacei, [xt, statet(:,2)]);
st   = Phi*cst(:,1) > Phi*cst(:,2) & xt + inj - kappau > sminu(1);   % switching decision of traditional sector producer     
xt(st==1) = xt(st==1) + inj(st==1) - kappau;             

% Sort savings into bins

% First, traditional

xtprime     = [xt(st==0); smint(1)];    % add entrants
[~, xtpbin] = histc(xtprime, edgest);   % bin by asset
xtpbin      = [xtpbin; (1:1:N)'];       % include all

% Second, unproductive

xuprime     = [xu(su==0); xt(st==1)];   % add entrants
[~, xupbin] = histc(xuprime, edgesu);   % bin by asset
xupbin  = [xupbin; (1:1:N)'];           % include all

% Third, productive

xpprime     = [xp; xu(su==1)];          % add entrants
[~, xppbin] = histc(xpprime, edgesp);   % bin by asset
xppbin  = [xppbin; (1:1:N)'];           % include all

% Iterate over measures

for t = 1 : T
    
    % first, traditional
    
ntold = nt;
nuold = nu;
npold = np;

ntprime     = [1/mu*ntold(st==0);    1-1/mu];
Ptprime     = [P(statet(st==0,2),:);  Perg'];
Ptprime     = bsxfun(@times, ntprime, Ptprime); 

Ptprime = [Ptprime; zeros(N,k)];    % make sure all guys included

Ptnew = zeros(N, k); 

for i = 1:k
    
   Ptnew(:,i) = accumarray(xtpbin,Ptprime(:,i));    % add all guys in a given bin
    
end

nt = reshape(Ptnew, N*k, 1); 

    % second, unproductive

nuprime = [1/mu*nuold(su==0); 1/mu*ntold(st==1)];
Puprime = [P(stateu(su==0,2),:); P(statet(st==1,2),:)];
Puprime = bsxfun(@times, nuprime, Puprime); 

Puprime = [Puprime; zeros(N,k)];    % make sure all guys included

Punew = zeros(N, k); 

for i = 1:k
    
   Punew(:,i) = accumarray(xupbin,Puprime(:,i));    % add all guys in a given bin
    
end

nu = reshape(Punew, N*k, 1); 


     % third, productive

npprime = [1/mu*npold; 1/mu*nuold(su==1)];
Ppprime = [P(statep(:,2),:); P(stateu(su==1,2),:)];
Ppprime = bsxfun(@times, npprime, Ppprime); 

Ppprime = [Ppprime; zeros(N,k)];    % make sure all guys included

Ppnew = zeros(N, k); 

for i = 1:k
    
   Ppnew(:,i) = accumarray(xppbin,Ppprime(:,i));    % add all guys in a given bin
    
end

np = reshape(Ppnew, N*k, 1); 

%fprintf('%4i %6.2e %6.2e %6.2e \n',[t, norm(nt-ntold), norm(nu-nuold), norm(np-npold)]);

if norm(nt-ntold) < 1e-10 && norm(nu-nuold) < 1e-10 && norm(np-npold) < 1e-10, break, end

end

% Productive Sector

rr   = funeval(crp,fspacep,statep);
Lp   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(egrid(statep(:,2))+phip);
Kp   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(egrid(statep(:,2))+phip);
Yp   = exp(egrid(statep(:,2))+phip).^(1-eta).*(Lp.^alpha.*Kp.^(1-alpha)).^eta;
Dp   = Yp - W.*Lp - (r+delta).*Kp - W*F;

Lt   = (eta/W)^(1/(1-eta))*exp(egrid(statep(:,2)));
Yt   = (eta/W)^(eta/(1-eta))*exp(egrid(statep(:,2)));
Dt   = pait*exp(egrid(statep(:,2)));

condmp = Dp >= Dt;

Kp = Kp.*condmp;
Lp = Lp.*condmp + Lt.*(~condmp);
Yp = Yp.*condmp + Yt.*(~condmp);
Dp = Dp.*condmp + Dt.*(~condmp);

Cp   = (1-theta*xai)*Dp + (1+r)*statep(:,1) - xp; 

% Unproductive Sector

rr   = funeval(cru,fspaceu,stateu);
Lu   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Ku   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Yu   = exp(egrid(stateu(:,2))+phiu).^(1-eta).*(Lu.^alpha.*Ku.^(1-alpha)).^eta;
Du   = Yu - W.*Lu - (r+delta).*Ku - W*F;  

Lt   = (eta/W)^(1/(1-eta))*exp(egrid(stateu(:,2)));
Yt   = (eta/W)^(eta/(1-eta))*exp(egrid(stateu(:,2)));
Dt   = pait*exp(egrid(stateu(:,2)));

condmu = Du >= Dt;

Ku = Ku.*condmu;
Lu = Lu.*condmu + Lt.*(~condmu);
Yu = Yu.*condmu + Yt.*(~condmu);
Du = Du.*condmu + Dt.*(~condmu);

Cu   = (1-theta*xai)*Du + (1+r)*stateu(:,1) - (xu + kappap*(su==1)); 

% Traditional Sector

Lt  = (eta/W)^(1/(1-eta))*exp(egrid(statet(:,2)));
Yt  = (eta/W)^(eta/(1-eta))*exp(egrid(statet(:,2)));
Dt  = Yt - W*Lt; 
Ct  = Dt + (1+r)*statet(:,1) - (xt + kappau*(st==1) - inj.*(st==1));   % we subtracted kappa from these guys, so add it back


% Substract equity from worker's savings

pu   = funeval(cu(:,3), fspaceu, stateu);
pp   = funeval(cu(:,3), fspacep, statep);

Eu   = 1/(1+r)*theta*xai*(Du + pu);
Ep   = 1/(1+r)*theta*xai*(Dp + pp);

Eq   = np'*Ep + nu'*Eu;

Ypa = np'*(Yp.*condmp);  % aggregate per capita output productive
Kpa = np'*(Kp.*condmp); 
Lpa = np'*(Lp.*condmp);

Yua = nu'*(Yu.*condmu);  % aggregate per capita output productive
Kua = nu'*(Ku.*condmu);  
Lua = nu'*(Lu.*condmu);

Lfa = np'*(F*condmp) + nu'*(F*condmu);

Yta = nt'*Yt + nu'*(Yu.*(~condmu)) + np'*(Yp.*(~condmp));  % traditional
Lta = nt'*Lt + nu'*(Lu.*(~condmu)) + np'*(Lp.*(~condmp));

Apa = np'*statep(:,1);  % assets in productive sector
Aua = nu'*stateu(:,1);  % assets in unproductive sector
Ata = nt'*statet(:,1);  % assets in traditional sector


Yma = Ypa + Yua;        % entire manufacturing
Lma = Lpa + Lua; 
Kma = Kpa + Kua; 

TFPm = Yma./(Lma.^alpha.*Kma.^(1-alpha)).^(eta);
TFPt = Yta./Lta.^eta;

TFPmbest = (np'*(exp(egrid(statep(:,2))+phip).*condmp) + nu'*(exp(egrid(stateu(:,2))+phiu).*condmu)).^(1-eta);



disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Kma,        Apa + Aua + Ata + Aw - Eq]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [Lta + Lma + Lfa,  Lbar]);
fprintf('\n');



D    = Kma - (Apa + Aua + Ata + Aw - Eq);    % Net Foreign Debt Position
Ya   = Yma + Yta; 
Ca   = Ya + (1 - delta - mu)*Kma - sum(np)*(mu-1)*kappap - sum(nu)*(mu-1)*kappau  - (1 + r - mu)*D; 
Ca2  = nt'*Ct + np'*Cp + nu'*Cu + W*Lbar + (1 + r - mu)*Aw;   % add consumption of all agents

fprintf('\n');
disp('Model Implications')
fprintf('\n');
fprintf('Fraction modern            = %9.3f \n',                sum(nu) + sum(np));
fprintf('Fraction modern & operate  = %9.3f \n',                sum(nu.*condmu) + sum(np.*condmp));

fprintf('Fraction output in modern  = %9.3f \n',          Yma/(Yma + Yta));
fprintf('Fraction labor in modern   = %9.3f \n',          Lma/(Lma + Lta));
fprintf('\n');

finvest = sum(np)*(mu-1)*kappap + sum(nu)*(mu-1)*kappau;

fprintf('Output                        = %9.3f \n',          Ya/(Lta+Lma+Lfa));   % per worker in case L doesnt clear
fprintf('Consumption                   = %9.3f \n',          Ca/(Lta+Lma+Lfa));
fprintf('Investment                    = %9.3f \n',          ((delta+mu-1)*Kma + finvest)/(Lta+Lma+Lfa));
fprintf('K to Y modern                 = %9.3f \n',          Kma/Yma);

fprintf('Fraction Fixed Investment     = %9.3f \n',          finvest/ ((delta+mu-1)*Kma + finvest));
fprintf('TFP Traditional               = %9.3f \n',          TFPt);
fprintf('TFP Modern                    = %9.3f \n',          TFPm);
fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPmbest/TFPm));
