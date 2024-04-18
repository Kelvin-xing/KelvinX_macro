Tinit = 100;       % periods to initialize distribution of assets
T     = 6;         % periods to simulate allocations
Nf    = 50000;    % # producers

randn('state', 100);
rand('state',  100);

zsim     = zeros(Nf,T);       % productivity
typesim  = zeros(Nf,T);       % end of period type: 1 = modern, 0 = traditional 
typeoldsim = zeros(Nf,T);     % beginning of period type
agesim   = zeros(Nf,T);       % age
agemsim  = zeros(Nf,T);       % age in modern sector
asim     = zeros(Nf,T);       % beginning of period assets
xsim     = zeros(Nf,T);       % savings for next period 

unif  = nodeunif(Nf,eps^(1/2),1-eps^(1/2));
unif  = unif(randperm(Nf));

jump  = nodeunif(Nf, eps, 1-eps) > 1/mu;
jump  = jump(randperm(Nf));

loset = nodeunif(Nf, eps, 1-eps) > rhot;   % get new productivity when transit from tradit. to modern
loset = loset(randperm(Nf)); 

losem = nodeunif(Nf, eps, 1-eps) > rhom;   % lose old productivity (and kappa) in modern sector
losem = losem(randperm(Nf)); 

% run a bit without saving data to get initial distribution

type  = zeros(Nf,1);    % start everyone in traditional sector

% Draw from unconditional distribution of z

wcum   = [0; cumsum(wnode)];
unif   = unif(randperm(Nf));
zdraw  = (bsxfun(@lt, unif, wcum(2:end)').*bsxfun(@ge, unif, wcum(1:end-1)'))*znode;
  
for t = 1 : Tinit+T

    if t==1
         
        z      = zdraw(randperm(Nf));
       
        age    = ones(Nf,1);
        agem   = zeros(Nf,1);
        
        state  = [smin(1)*ones(Nf,1), z]; 
        a      = state(:,1);
        x      = funeval(cx(:,2), fspace, state);

        Phi    = funbas(fspace, [x,z]);        % Evaluate value of switching/staying next period give current x
        
        type       = Phi*cst(:,1)>Phi*cst(:,2) & x - kappa*exp(state(:,2)) > smin(1);
        x(type==1) = x(type==1) - kappa*exp(state(type==1,2)); 
        
        sw = (type==1); % keep track of who switched
        
    else

        jump   = jump(randperm(Nf));
        loset  = loset(randperm(Nf));
        losem  = losem(randperm(Nf)); 
        
        zdraw  = zdraw(randperm(Nf)); 
        
        z(sw==1)    = (1-loset(sw==1)).*z(sw==1) + loset(sw==1).*zdraw(sw==1);
        z(sw~=1)    = (1-losem(sw~=1)).*z(sw~=1) + losem(sw~=1).*zdraw(sw~=1);
        type        = type.*(~losem); 
        
        type   = type.*(~jump);
        age    = age.*(~jump) + 1;
        agem   = (agem + 1).*(type>=1);
        
        x      = x.*(~jump) + smin(1)*(jump);
        z      = z.*(~jump) + zdraw.*(jump);
        
      typeold  = type;   % save to determine who pays the fixed cost
        
      sw  = zeros(Nf,1); % reset switching indicator
      
   % Compute savings rules for both types
  
        if any(typeold == 1)
            
            state      = [x(typeold == 1), z(typeold == 1)];
            a(typeold == 1) = state(:,1);
            x(typeold == 1) = max(funeval(cx(:,1), fspace, state), smin(1)); 

            type(typeold == 1) = type(typeold == 1);

        end
        
        
        if any(typeold == 0)
                        
            state         = [x(typeold == 0), z(typeold==0)];
            a(typeold==0) = state(:,1);
            x(typeold==0) = max(funeval(cx(:,2), fspace, state), smin(1)); 
           
            Phi           = funbas(fspace, [x(typeold==0), z(typeold==0)]);
            
            type(typeold==0) = Phi*cst(:,1)>Phi*cst(:,2) & x(typeold==0) - kappa*exp(z(typeold==0)) > smin(1);
            x(type == 1 & typeold ==0) = x(type==1 & typeold ==0) - kappa*exp(z(type==1 & typeold==0));

            sw(typeold==0) = type(typeold==0); 
            
        end
              
          if t > Tinit
             
              typesim(:,t-Tinit) = type;
              agesim(:,t-Tinit)  = age;
              agemsim(:,t-Tinit) = agem;
              asim(:,t-Tinit)    = a;
              zsim(:,t-Tinit)    = z;
              xsim(:,t-Tinit)    = x + kappa.*exp(z).*(type == 1 & typeold == 0);    % convert back to savings prior to paying fixed cost
              typeoldsim(:,t-Tinit) = typeold;

          end
    end
end

sw   = (typesim' == 1 & typeoldsim' == 0);
a    = asim';
x    = xsim';
age  = agesim';
agem = agemsim';
z    = zsim';
type = typeoldsim';     % beginning of period type   
rm   = reshape(funeval(cr,fspace,[a(:), z(:)]), size(a,1), size(a,2));
rr   = rm.*(type==1) + r.*(type==0);

Lm   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(z + phi);
Km   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(z + phi);
Ym   = exp(z+phi).^(1-eta).*(Lm.^alpha.*Km.^(1-alpha)).^eta;
Dm   = Ym - W.*Lm - (r+delta).*Km;  

Lt   = (eta/W)^(1/(1-eta))*exp(z);
Yt   = (eta/W)^(eta/(1-eta))*exp(z);
Dt   = Yt - W.*Lt;

Y    = Ym.*(type==1) + Yt.*(type==0);
K    = Km.*(type==1);
L    = Lm.*(type==1) + Lt.*(type==0);
D    = Dm.*(type==1) + Dt.*(type==0);

cons = D + (1+r)*a - x;                    % consumption, (recall x is total savings)

% Frictionless allocations

Lmf   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(z + phi);
Kmf   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(z + phi);
Ymf   = exp(z+phi).^(1-eta).*(Lmf.^alpha.*Kmf.^(1-alpha)).^eta;

Ltf   = (eta/W)^(1/(1-eta))*exp(z);
Ytf   = (eta/W)^(eta/(1-eta))*exp(z);

Yf    = Ymf.*(type==1) + Ytf.*(type==0);
Kf    = Kmf.*(type==1);
Lf    = Lmf.*(type==1) + Ltf.*(type==0);

A = a;

Debt   = K - A;

% three ways to compute TFP losses in Modern Sector (TFP calculation only
% up to a scalar of normalization, may differ from number in ergodic.m)

zt = z + phi.*(type==1);   % total efficiency (exogenous and endogenous)

TFP      = mean(Y(type>0))./(mean(L(type>0)).^alpha*mean(K(type>0)).^(1-alpha))^eta;
TFPbest  = mean(Yf(type>0))./(mean(Lf(type>0)).^alpha*mean(Kf(type>0)).^(1-alpha))^eta;
TFP2best = mean(exp(zt(type>0))).^(1-eta);

wedgel = (rr+delta).^(-(1-alpha)*eta/(1-eta)).*mean(exp(zt(type>0)))./mean((rr(type>0)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(zt(type>0)));
wedgek = (rr+delta).^((alpha*eta-1)/(1-eta)).*mean(exp(zt(type>0)))./mean((rr(type>0)+delta).^((alpha*eta-1)/(1-eta)).*exp(zt(type>0)));

TFPratio = mean(wedgel(type>0).^(alpha*eta).*wedgek(type>0).^((1-alpha)*eta).*exp(zt(type>0)))./mean(exp(zt(type>0)));

agem(1,:) = agem(1,:).*(type(1,:)==1);

for t = 2:size(agem,1)
    
    agem(t,:) = (agem(t-1,:) + 1).*(type(t,:)==1);
    
end


% autocorrelation Y for firms in modern sector

Yc   = vec(Y(2:end,:));
Yp   = vec(Y(1:end-1,:));
Lc   = vec(L(2:end,:));
Lp   = vec(L(1:end-1,:));
Kc   = vec(K(2:end,:));
Kp   = vec(K(1:end-1,:));

flag = vec(agem(2:end,:)>1); 

acorY = corr([log(Yc(flag)),log(Yp(flag))]);
acorL = corr([log(Lc(flag)),log(Lp(flag))]);
acorK = corr([log(Kc(flag)),log(Kp(flag))]);

Yc = vec(Y(3:end,:));
Yp = vec(Y(1:end-2,:));
Lc = vec(L(3:end,:));
Lp = vec(L(1:end-2,:));
Kc = vec(K(3:end,:));
Kp = vec(K(1:end-2,:));

flag = vec(agem(3:end,:)>2); 

acorY2 = corr([log(Yc(flag)), log(Yp(flag))]);
acorL2 = corr([log(Lc(flag)), log(Lp(flag))]);
acorK2 = corr([log(Kc(flag)), log(Kp(flag))]);

Yc = vec(Y(4:end,:));
Yp = vec(Y(1:end-3,:));
Lc = vec(L(4:end,:));
Lp = vec(L(1:end-3,:));
Kc = vec(K(4:end,:));
Kp = vec(K(1:end-3,:));

flag = vec(agem(4:end,:)>3); 

acorY3 = corr([log(Yc(flag)), log(Yp(flag))]);
acorL3 = corr([log(Lc(flag)), log(Lp(flag))]);
acorK3 = corr([log(Kc(flag)), log(Kp(flag))]);

Yc = vec(Y(5:end,:));
Yp = vec(Y(1:end-4,:));
Lc = vec(L(5:end,:));
Lp = vec(L(1:end-4,:));
Kc = vec(K(5:end,:));
Kp = vec(K(1:end-4,:));

flag = vec(agem(5:end,:)>4); 

acorY4 = corr([log(Yc(flag)), log(Yp(flag))]);
acorL4 = corr([log(Lc(flag)), log(Lp(flag))]);
acorK4 = corr([log(Kc(flag)), log(Kp(flag))]);

Yc = vec(Y(6:end,:));
Yp = vec(Y(1:end-5,:));
Lc = vec(L(6:end,:));
Lp = vec(L(1:end-5,:));
Kc = vec(K(6:end,:));
Kp = vec(K(1:end-5,:));

flag = vec(agem(6:end,:)>5); 

acorY5 = corr([log(Yc(flag)), log(Yp(flag))]);
acorL5 = corr([log(Lc(flag)), log(Lp(flag))]);
acorK5 = corr([log(Kc(flag)), log(Kp(flag))]);

% growth rates:

Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));
Lc = vec(L(2:end,:));
Lp = vec(L(1:end-1,:));
Kc = vec(K(2:end,:));
Kp = vec(K(1:end-1,:));

flag = vec(agem(2:end,:)>1); 

dY  = vec(log(Yc(flag)./Yp(flag))); 
dL  = vec(log(Lc(flag)./Lp(flag)));
dK  = vec(log(Kc(flag)./Kp(flag)));

% firm-level 'finance' statistics

YK3 = mean(log(Y(agem>10&type>0)./(K(agem>10&type>0))));
YK1 = mean(log(Y(agem<=5&type>0)./(K(agem<=5&type>0))))-YK3;
YK2 = mean(log(Y(agem>5&agem<=10&type>0)./(K(agem>5&agem<=10&type>0))))-YK3;

YKa3 = log(sum(Y(agem>10&type>0))./sum(K(agem>10&type>0)));
YKa1 = log(sum(Y(agem<=5&type>0))./sum(K(agem<=5&type>0))) - YKa3;
YKa2 = log(sum(Y(agem>5&agem<=10&type>0))./sum(K(agem>5&agem<=10&type>0))) - YKa3;

LK3 = mean(log(L(agem>10&type>0)./(K(agem>10&type>0))));
LK1 = mean(log(L(agem<=5&type>0)./(K(agem<=5&type>0))))-LK3;
LK2 = mean(log(L(agem>5&agem<=10&type>0)./(K(agem>5&agem<=10&type>0))))-LK3;

LKa3 = log(sum(L(agem>10&type>0))./sum(K(agem>10&type>0)));
LKa1 = log(sum(L(agem<=5&type>0))./sum(K(agem<=5&type>0))) - LKa3;
LKa2 = log(sum(L(agem>5&agem<=10&type>0))./sum(K(agem>5&agem<=10&type>0))) - LKa3;

Yc    = vec(Y(2:end,:));
Yp    = vec(Y(1:end-1,:));
Lc    = vec(L(2:end,:));
Lp    = vec(L(1:end-1,:));
Kc    = vec(K(2:end,:));
Kp    = vec(K(1:end-1,:));

flag1 = vec(agem(2:end,:)>1 & agem(2:end,:)<=5); 
flag2 = vec(agem(2:end,:)>5 & agem(2:end,:)<=10 );
flag3 = vec(agem(2:end,:)>10);

dY1   = mean(log(Yc(flag1)./Yp(flag1))); 
dY2   = mean(log(Yc(flag2)./Yp(flag2))); 
dY3   = mean(log(Yc(flag3)./Yp(flag3))); 

dL1   = mean(log(Lc(flag1)./Lp(flag1))); 
dL2   = mean(log(Lc(flag2)./Lp(flag2))); 
dL3   = mean(log(Lc(flag3)./Lp(flag3))); 

dK1   = mean(log(Kc(flag1)./Kp(flag1))); 
dK2   = mean(log(Kc(flag2)./Kp(flag2))); 
dK3   = mean(log(Kc(flag3)./Kp(flag3))); 

dYa1   = log(sum(Yc(flag1))/sum(Yp(flag1))); 
dYa2   = log(sum(Yc(flag2))/sum(Yp(flag2))); 
dYa3   = log(sum(Yc(flag3))/sum(Yp(flag3)));  

dLa1   = log(sum(Lc(flag1))/sum(Lp(flag1))); 
dLa2   = log(sum(Lc(flag2))/sum(Lp(flag2))); 
dLa3   = log(sum(Lc(flag3))/sum(Lp(flag3))); 

dKa1   = log(sum(Kc(flag1))/sum(Kp(flag1))); 
dKa2   = log(sum(Kc(flag2))/sum(Kp(flag2))); 
dKa3   = log(sum(Kc(flag3))/sum(Kp(flag3))); 

% Pin down variance of z needed to match the size distributions

ac1 = acorY; 
ac3 = acorY3;
ac5 = acorY5; 
sY = std(log(Y(type>0)));

fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

data = sort(Ym(:),1);

f10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);

disp('Moments');

fprintf('\n');
disp('                   1st col: model,      2nd col: data');
fprintf('\n');
fprintf('s.d.  dlog(Y)                 = %9.2f  %9.2f \n',       [std(dY(:)),           sdYdata]);
fprintf('s.d.   log(Y)                 = %9.2f  %9.2f \n',       [sY,                   sYdata]);
fprintf('Fraction revenue top 10 perc  = %9.2f  %9.2f \n',       [f10, 0.84]);
fprintf('Modern to Trad L ratio        = %9.2f  %9.2f \n',       [mean(L(type==1))/mean(L(type==0)), 5])
fprintf('\n');
fprintf('cor y y1                      = %9.2f  %9.2f \n',       [ac1,                  ac1Ydata]);
fprintf('cor y y3                      = %9.2f  %9.2f \n',       [ac3,                  ac3Ydata]);
fprintf('cor y y5                      = %9.2f  %9.2f \n',       [ac5,                  ac5Ydata]);

fprintf('\n');
fprintf('aggregate Debt to GDP         = %9.2f  %9.2f \n',       [sum(Debt(type>0))./sum(Y(type>0)),       DtoY]); 
fprintf('aggregate posit. Debt to GDP  = %9.2f  %9.2f \n',       [sum(Debt(Debt>0 & type>0))./sum(Y(type>0)),       DtoY]); 
fprintf('\n');
finvest = sum(sw(:).*kappa.*exp(z(:)));

fprintf('Fixed Inv  to total Y, perc   = %9.1f  %9.1f \n',       [finvest/sum(Y(agem>=1))*100, 4.6]);
fprintf('\n');
fprintf('fraction producers 1-10       = %9.2f  %9.2f \n',       [sum(agem(:)>=1 & agem(:)< 11)/sum(agem(:)>=1),             0.51 + .26]);
fprintf('\n');

%fprintf('Variance exog component       = %9.2f \n',       varz);
 if 1

fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('Jiwoon, ignore rest')
fprintf('\n');
fprintf('\n');
fprintf('\n');


disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(K(:)),   mean(A(:))]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(L(:)),   Lbar]);
fprintf('\n');

fprintf('\n');
disp('Model Implications')
fprintf('\n');
fprintf('Fraction modern               = %9.3f \n',          mean(type(:)==1));
fprintf('Fraction output in modern     = %9.3f \n',          sum(Y(type>0))/sum(Y(:)));
fprintf('Fraction labor in modern      = %9.3f \n',          sum(L(type>0))/sum(L(:)));
fprintf('Fraction new entrants         = %9.3f \n',          sum(agem(:)==1)/sum(agem(:)>=1));

fprintf('\n');

Da = mean(K(:)) - mean(A(:));       % Aggregate debt position
finvest = mean(sw(:).*kappa.*exp(z(:)));

fprintf('Output                        = %9.3f \n',          mean(Y(:))/mean(L(:)));
fprintf('Consumption                   = %9.3f \n',          (mean(Y(:)) - (delta+mu-1)*mean(K(:)) - (1+r-mu)*Da  - finvest)/mean(L(:)));
fprintf('Investment                    = %9.3f \n',          ((delta + mu - 1)*mean(K(:)) + finvest)/mean(L(:)));
fprintf('Fraction Fixed Investment     = %9.3f \n',          finvest/ ((delta + mu - 1)*mean(K(:)) + finvest));
fprintf('TFP Traditional               = %9.3f \n',          sum(Y(type==0))/sum(L(type==0)).^(eta)*numel(Y)^(eta-1));   % express in per capita units, comparable to ergodic
fprintf('TFP Modern                    = %9.3f \n',          sum(Y(type>0))./(sum(L(type>0)).^alpha*sum(K(type>0)).^(1-alpha))^eta*numel(Y)^(eta-1)); % express in per capita units, comparable to ergodic
fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPbest/TFP)*100);


fprintf('\n');


ract = rr(type>0);

fprintf('\n');
fprintf('Fraction constrained         = %10.2f \n', mean(   ract > r + 0.0025));
fprintf('Median r if constrained      = %10.2f \n', median( ract(ract > r + 0.0025)-r));
fprintf('Iqr r if constrained         = %10.2f \n', iqr(    ract(ract > r + 0.0025)-r));
fprintf('90 pctile r if constrained   = %10.2f \n', prctile(ract(ract > r + 0.0025)-r,90));
fprintf('99 pctile r if constrained   = %10.2f \n', prctile(ract(ract > r + 0.0025)-r,99));
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

disp('Other moments')

fprintf('\n');
fprintf('fraction producers 1-5        = %9.2f  %9.2f \n',       [sum(agem(:)>=1 & agem(:)< 6)/sum(agem(:)>=1),             0.51]);
fprintf('fraction producers 6-10       = %9.2f  %9.2f \n',       [sum(agem(:)>=6 & agem(:)< 11)/sum(agem(:)>=1),            0.26]);

fprintf('fraction output 1-5           = %9.2f  %9.2f \n',       [sum(Y(agem(:)>=1 & agem(:)< 6))/sum(Y(agem(:)>=1)),             0.20]);
fprintf('fraction output 6-10          = %9.2f  %9.2f \n',       [sum(Y(agem(:)>=6 & agem(:)< 11))/sum(Y(agem(:)>=1)),            0.20]);

fprintf('fraction capital 1-5           = %9.2f  %9.2f \n',       [sum(K(agem(:)>=1 & agem(:)< 6))/sum(K(agem(:)>=1)),             0.20]);
fprintf('fraction capital 6-10          = %9.2f  %9.2f \n',       [sum(K(agem(:)>=6 & agem(:)< 11))/sum(K(agem(:)>=1)),            0.20]);

fprintf('\n');


fprintf('aggr DlogY 1-5                = %9.2f  %9.2f \n',     [dYa1,     0.199]);
fprintf('aggr DlogY 6-10               = %9.2f  %9.2f \n',     [dYa2,     0.089]);
fprintf('aggr DlogL 1-5                = %9.2f  %9.2f \n',     [dLa1,     0.162]);
fprintf('aggr DlogL 6-10               = %9.2f  %9.2f \n',     [dLa2,     0.068]);
fprintf('aggr DlogK 1-5                = %9.2f  %9.2f \n',     [dKa1,     0.249]);
fprintf('aggr DlogK 6-10               = %9.2f  %9.2f \n',     [dKa2,     0.040]);
fprintf('\n');
fprintf('mean  YK 1-5                  = %9.2f  %9.2f \n',     [YK1,      0.09]);
fprintf('mean  YK 6-10                 = %9.2f  %9.2f \n',     [YK2,      0.03]);
fprintf('mean  LK 1-5                  = %9.2f  %9.2f \n',     [LK1,      0.11]);
fprintf('mean  LK 6-10                 = %9.2f  %9.2f \n',     [LK2,      0.03]);
fprintf('\n');
fprintf('aggr  YK 1-5                  = %9.2f  %9.2f \n',     [YKa1,    -0.03]);
fprintf('aggr  YK 6-10                 = %9.2f  %9.2f \n',     [YKa2,     0.03]);
fprintf('aggr  LK 1-5                  = %9.2f  %9.2f \n',     [LKa1,    -0.01]);
fprintf('aggr  LK 6-10                 = %9.2f  %9.2f \n',     [LKa2,     0.03]);
fprintf('\n');

varL  = var(log(L(type>0)));
sdL   = sqrt(varL);
acL1  = acorL;
acL3  = acorL3;
acL5  = acorL5;

varK  = var(log(K(type>0)));
sdK   = sqrt(varK);
acK1  = acorK;
acK3  = acorK3;
acK5  = acorK5;

fprintf('\n');
fprintf('s.d.  dlog(L)                 = %9.2f  %9.2f \n',     [std(dL),           0.49]);
fprintf('s.d.   log(L)                 = %9.2f  %9.2f \n',     [sdL,               1.22]);
fprintf('\n');
fprintf('s.d.  dlog(K)                 = %9.2f  %9.2f \n',     [std(dK),           0.57]);
fprintf('s.d.   log(K)                 = %9.2f  %9.2f \n',     [sdK,               1.44]);
fprintf('\n');
fprintf('cor L L1                      = %9.2f  %9.2f \n',     [acL1,              0.92]);
fprintf('cor L L3                      = %9.2f  %9.2f \n',     [acL3,              0.89]);
fprintf('cor L L5                      = %9.2f  %9.2f \n',     [acL5,              0.86]);
fprintf('\n');
fprintf('cor K K1                      = %9.2f  %9.2f \n',     [acK1,              0.92]);
fprintf('cor K K3                      = %9.2f  %9.2f \n',     [acK3,              0.89]);
fprintf('cor K K5                      = %9.2f  %9.2f \n',     [acK5,              0.85]);
fprintf('\n');

fprintf('\n');

 end