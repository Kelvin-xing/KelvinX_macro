T = 1000;  % periods we wait for weights to converge

N = 5001;

agridu = nodeunif(N, sminu(1), smaxu(1));
agridt = nodeunif(N, smint(1), smaxt(1));

edgesu = [agridu(1); (agridu(2:end-1) + agridu(1:end-2))/2; agridu(end)];
edgest = [agridt(1); (agridt(2:end-1) + agridt(1:end-2))/2; agridt(end)];

stateu = gridmake(agridu, (1:1:k)');
statet = gridmake(agridt, (1:1:k)');

Nu = size(stateu, 1); 
Nt = size(statet, 1);

nu = 0.70*1/Nu*ones(Nu,1);
nt = 0.30*1/Nt*ones(Nt,1);

% Decision rules

xu   = funeval(cxu, fspaceu, stateu); xu = min(max(xu, sminu(1)),smaxu(1)); % don't allow outside bounds 
xt   = funeval(cxt, fspacet, statet); xt = min(max(xt, smint(1)),smaxt(1)); % don't allow outside bounds 


Phi  = funbas(fspacet,        [xt, statet(:,2)]);
st   = Phi*cst(:,1) > Phi*cst(:,2) & xt  - kappau > sminu(1);   % switching decision of traditional sector producer     
xt(st==1) = xt(st==1) - kappau;             

% Sort savings into bins

% First, traditional

xtprime     = [xt(st==0); smint(1)];    % add entrants
[~, xtpbin] = histc(xtprime, edgest);   % bin by asset
xtpbin      = [xtpbin; (1:1:N)'];       % include all

% Second, unproductive

xuprime     = [xu; xt(st==1)];          % add entrants
[~, xupbin] = histc(xuprime, edgesu);   % bin by asset
xupbin  = [xupbin; (1:1:N)'];           % include all


% Iterate over measures

for t = 1 : T
    
    % first, traditional
    
ntold = nt;
nuold = nu;

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

nuprime = [1/mu*nuold; 1/mu*ntold(st==1)];
Puprime = [P(stateu(:,2),:); P(statet(st==1,2),:)];
Puprime = bsxfun(@times, nuprime, Puprime); 

Puprime = [Puprime; zeros(N,k)];    % make sure all guys included

Punew = zeros(N, k); 

for i = 1:k
    
   Punew(:,i) = accumarray(xupbin,Puprime(:,i));    % add all guys in a given bin
    
end

nu = reshape(Punew, N*k, 1); 

%fprintf('%4i %6.2e %6.2e  \n',[t, norm(nt-ntold), norm(nu-nuold)]);

if norm(nt-ntold) < 1e-10 && norm(nu-nuold) < 1e-10, break, end

end

% Unproductive Sector

rr   = funeval(cru,fspaceu,stateu);
Lu   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Ku   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Yu   = exp(egrid(stateu(:,2))+phiu).^(1-eta).*(Lu.^alpha.*Ku.^(1-alpha)).^eta;
Du   = Yu - W.*Lu - (r+delta).*Ku;  
Cu   = Du + (1+r)*stateu(:,1) - xu; 

% Traditional Sector

Lt  = (eta/W)^(1/(1-eta))*exp(egrid(statet(:,2)));
Yt  = (eta/W)^(eta/(1-eta))*exp(egrid(statet(:,2)));
Dt  = Yt - W*Lt; 
Ct  = Dt + (1+r)*statet(:,1) - (xt + kappau*(st==1));   % we subtracted kappa from these guys, so add it back

Yua = nu'*Yu;  % aggregate per capita output productive
Kua = nu'*Ku;  
Lua = nu'*Lu;

Yta = nt'*Yt;  % traditional
Lta = nt'*Lt;

Aua = nu'*stateu(:,1);  % assets in unproductive sector
Ata = nt'*statet(:,1);  % assets in traditional sector

TFPu = Yua./(Lua.^alpha.*Kua.^(1-alpha)).^(eta);
TFPt = Yta./Lta.^eta;

TFPubest = (nu'*exp(egrid(stateu(:,2))+phiu)).^(1-eta);

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Kua,        Aua + Ata]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [Lta + Lua,  Lbar]);
fprintf('\n');

D    = Kua - (Aua + Ata);    % Net Foreign Debt Position
Ya   = Yua + Yta; 
Ca   = Ya + (1 - delta - mu)*Kua - sum(nu)*(mu-1)*kappau  - (1 + r - mu)*D; 
Ca2  = nt'*Ct + nu'*Cu + W*Lbar;   % add consumption of all agents

fprintf('\n');
disp('Model Implications')
fprintf('\n');
fprintf('Fraction modern               = %9.3f \n',          sum(nu));
fprintf('Fraction output in modern     = %9.3f \n',          Yua/(Yua + Yta));
fprintf('Fraction labor in modern      = %9.3f \n',          Lua/(Lua + Lta));
fprintf('\n');

finvest = + sum(nu)*(mu-1)*kappau;

fprintf('Output                        = %9.3f \n',          Ya/(Lta+Lua));   % per worker in case L doesnt clear
fprintf('Consumption                   = %9.3f \n',          Ca/(Lta+Lua));
fprintf('Investment                    = %9.3f \n',          ((delta+mu-1)*Kua + finvest)/(Lta+Lua));
fprintf('Fraction Fixed Investment     = %9.3f \n',          finvest/ ((delta+mu-1)*Kua + finvest));
fprintf('TFP Traditional               = %9.3f \n',          TFPt);
fprintf('TFP Modern                    = %9.3f \n',          TFPu);
fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPubest/TFPu));


save steady_state N agridu agridt edgesu edgest Nu Nt nu nt stateu statet
xinit = [Yua; Kua; Lua; Yta; Lta; Aua; Ata; TFPu; TFPt; TFPubest; D; Ya; Ca];

save xinit xinit