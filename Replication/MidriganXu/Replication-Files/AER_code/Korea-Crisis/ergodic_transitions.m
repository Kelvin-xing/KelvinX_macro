load steady_state; % load N agridu agridt edgesu edgest Nu Nt nu nt stateu statet

ntsave = zeros(size(nt,1),  T);
nusave = zeros(size(nt,1),  T);

fprintf('\n')
disp('Computing transitions given decision rules')

for t = 1 : T
  
ntsave(:,t) = nt;
nusave(:,t) = nu;
    
    cxu = cxu_irf(:,:,t); 
    cxt = cxt_irf(:,:,t);
    cst = cst_irf(:,:,t);
    
    
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
    
% first, traditional
    
ntprime     = [1/mu*ntsave(st==0,t);    1-1/mu];
Ptprime     = [P(statet(st==0,2),:);     Perg'];
Ptprime     = bsxfun(@times, ntprime, Ptprime); 

Ptprime = [Ptprime; zeros(N,k)];    % make sure all guys included

Ptnew = zeros(N, k); 

for i = 1:k
    
   Ptnew(:,i) = accumarray(xtpbin, Ptprime(:,i));    % add all guys in a given bin
    
end

nt = reshape(Ptnew, N*k, 1); 

    % second, unproductive

nuprime = [1/mu*nusave(:,t); 1/mu*ntsave(st==1,t)];
Puprime = [P(stateu(:,2),:); P(statet(st==1,2),:)];
Puprime = bsxfun(@times, nuprime, Puprime); 

Puprime = [Puprime; zeros(N,k)];    % make sure all guys included

Punew = zeros(N, k); 

for i = 1:k
    
   Punew(:,i) = accumarray(xupbin,Puprime(:,i));    % add all guys in a given bin
    
end

nu = reshape(Punew, N*k, 1); 

end

Yua       = zeros(T,1);
Kua       = zeros(T,1);
Lua       = zeros(T,1); 
Yta       = zeros(T,1);
Lta       = zeros(T,1);
Aua       = zeros(T,1);
Ata       = zeros(T,1);
TFPu      = zeros(T,1);
TFPt      = zeros(T,1);
TFPubest  = zeros(T,1);
Ca        = zeros(T,1);
D         = zeros(T,1);       % Net Foreign Debt Position
Ya        = zeros(T,1);
num       = zeros(T,1);       % measure modern sector

for t = 1 : T

% Unproductive Sector

cru = cru_irf(:,:,t);
W   = Wsave(t);

rr   = funeval(cru,fspaceu,stateu);
Lu   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Ku   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(egrid(stateu(:,2))+phiu);
Yu   = exp(egrid(stateu(:,2))+phiu).^(1-eta).*(Lu.^alpha.*Ku.^(1-alpha)).^eta;
Du   = Yu - W.*Lu - (r+delta).*Ku;  

% Traditional Sector

Lt  = (eta/W)^(1/(1-eta))*exp(egrid(statet(:,2)));
Yt  = (eta/W)^(eta/(1-eta))*exp(egrid(statet(:,2)));
Dt  = Yt - W*Lt; 
    
Yua(t) = nusave(:,t)'*Yu;  % aggregate per capita output productive
Kua(t) = nusave(:,t)'*Ku;  
Lua(t) = nusave(:,t)'*Lu;

Yta(t) = ntsave(:,t)'*Yt;  % traditional
Lta(t) = ntsave(:,t)'*Lt;

Aua(t) = nusave(:,t)'*stateu(:,1);  % assets in unproductive sector
Ata(t) = ntsave(:,t)'*statet(:,1);  % assets in traditional sector

TFPu(t) = Yua(t)./(Lua(t).^alpha.*Kua(t).^(1-alpha)).^(eta);
TFPt(t) = Yta(t)./Lta(t).^eta;

TFPubest(t) = (nusave(:,t)'*exp(egrid(stateu(:,2))+phiu)).^(1-eta);

D(t)    = Kua(t) - (Aua(t) + Ata(t));    % Net Foreign Debt Position
Ya(t)   = Yua(t) + Yta(t); 
num(t) = sum(nusave(:,t));

end

Ynet = zeros(T,1);

for t = 1 : T-1
    
    Ca(t) = Ya(t) + (1-delta)*Kua(t) - sum(nusave(:,t))*(mu-1)*kappau - mu*Kua(t+1) + mu*D(t+1)  - (1 + r)*D(t);
    Ynet(t) = Ya(t) + (1-delta)*Kua(t) - sum(nusave(:,t))*(mu-1)*kappau - mu*Kua(t+1); 
end

Ca(T) = Ca(T-1);
Ynet(T) = Ynet(T-1);

disp(norm(Lta + Lua - 1))

Wsave = Wsave.*(Lua + Lta).^(1/5);

save Wsave Wsave;

load xinit; % Yua, Kua, Lua, Yta Lta Aua Ata TFPu TFPt TFPubest D Ya Ca

Yua  = [xinit(1); Yua];
Kua  = [xinit(2); Kua];
Lua  = [xinit(3); Lua];
Yta  = [xinit(4); Yta];
Lta  = [xinit(5); Lta];
Aua  = [xinit(6); Aua];
Ata  = [xinit(7); Ata];
TFPu = [xinit(8); TFPu];
TFPt = [xinit(9); TFPt];
TFPubest = [xinit(10); TFPubest];
D    = [xinit(11); D];
Ya   = [xinit(12); Ya];
Ca   = [xinit(13); Ca];
num  = [num(1); num];


close all
set(gcf,'DefaultLineLineWidth',2);
subplot(2,2,1)
hold on
plot((0:1:10)',[num(1:11), Lua(1:11)./(Lua(1:11) + Lta(1:11))]);
legend('measure producers modern', 'employment modern')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',18,'LineWidth',2)

subplot(2,2,2)
hold on
plot((0:1:10)', log([Yua(1:11)/Yua(1), Lua(1:11)/Lua(1), Kua(1:11)/Kua(1)]))
legend('output modern', 'employment modern', 'capital modern')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',18,'LineWidth',2)

subplot(2,2,3)
hold on
plot((0:1:10)', log([TFPu(1:11)/TFPu(1), TFPubest(1:11)/TFPubest(1)]))
legend('TFP modern', 'First Best TFP Modern')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',18,'LineWidth',2)

subplot(2,2,4)
hold on
plot((0:1:10)',log(TFPubest(1:11)./TFPu(1:11)))
legend('TFP loss misallocation')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',18,'LineWidth',2)
