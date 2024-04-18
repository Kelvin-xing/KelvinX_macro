Tinit = 50;       % periods to initialize distribution of assets
T     = 6;        % periods to simulate allocations
Nf    = 50000;    % # producers

randn('state', 100);
rand('state',  100);

eesim    = zeros(Nf,T);       % transitory component
typesim  = zeros(Nf,T);       % end of period type: 2 = productive, 1 = unproductive, 0 = traditional 
typeoldsim = zeros(Nf,T);     % beginning of period type
agesim   = zeros(Nf,T);       % age
agemsim  = zeros(Nf,T);       % age in modern sector
asim     = zeros(Nf,T);       % beginning of period assets
xsim     = zeros(Nf,T);        

Pcum     = P*triu(ones(size(P)));   % cumulative ergodic distributon

unif = nodeunif(Nf,eps^(1/2),1-eps^(1/2));
unif = unif(randperm(Nf));

ppi       = [0; cumsum(Perg)];
[n,bin]   = histc(unif,ppi);

esim      = zeros(Nf,T+Tinit);
esim(:,1) = bin;

for t = 2 : Tinit + T
   
  unif = unif(randperm(Nf));
  Pnew = P(esim(:,t-1),:);
  Pcum = [zeros(Nf,1), cumsum(Pnew,2)];
  esim(:,t) = ((repmat(unif,1,k)<Pcum(:,2:end)).*(repmat(unif,1,k)>Pcum(:,1:end-1)))*(1:1:k)';
  
end

jump  = nodeunif(Nf, eps, 1-eps) > 1/mu;
jump  = jump(randperm(Nf));

% run a bit without saving data to get initial distribution

type  = zeros(Nf,1);    % start everyone in traditional sector

for t = 1 : Tinit+T

    if t==1
         
        e      = esim(:,t);
        age    = ones(Nf,1);
        agem   = zeros(Nf,1);
        
        state  = [smint(1)*ones(Nf,1), e]; 
        a      = state(:,1);
        x      = funeval(cxt,fspacet,state);

        Phi    = funbas(fspacet,[x,e]);        % Evaluate value of switching/staying next period give current x
        
        type   = Phi*cst(:,1)>Phi*cst(:,2) & x - kappau > sminu(1);
        x(type==1) = x(type==1) - kappau; 
                
    else

        jump   = jump(randperm(Nf));
       
        e      = esim(:,t);                  % "exit" is random, so ok to use incumbent's e as draw from ergodic
        type   = type.*(~jump);

        age    = age.*(~jump) + 1;
        agem   = (agem + 1).*(type>=1);
        
        x      = x.*(~jump) + smint(1)*(jump);
 
  typeold    = type;   % save to determine who pays the fixed cost
        
        
        if any(typeold == 2)
            
           state           = [x(typeold == 2), e(typeold==2)];
           a(typeold == 2) = state(:,1);
           x(typeold == 2) = max(funeval(cxp, fspacep, state), sminp(1));
            
        end
        
        
        if any(typeold == 1)
            
            state      = [x(typeold == 1), e(typeold == 1)];
            a(typeold == 1) = state(:,1);
            x(typeold == 1) = max(funeval(cxu, fspaceu, state), sminu(1)); 
           
            Phi        = funbas(fspaceu, [x(typeold == 1), e(typeold == 1)]);

            type(typeold == 1) = type(typeold == 1) + (Phi*csu(:,1)>Phi*csu(:,2) & x(typeold == 1) - kappap > sminp(1));
            x(type == 2 & typeold ==1) = x(type==2 & typeold ==1) - kappap;

        end
        
        
        if any(typeold == 0)
                        
            state         = [x(typeold == 0), e(typeold==0)];
            a(typeold==0) = state(:,1);
            x(typeold==0) = max(funeval(cxt, fspacet, state), smint(1)); 
           
            Phi           = funbas(fspacet, [x(typeold==0), e(typeold==0)]);

            inj = funeval(cinj, fspacei, [x(typeold==0), e(typeold==0)]); 
            
            type(typeold==0) = Phi*cst(:,1)>Phi*cst(:,2) & x(typeold==0) + inj - kappau > sminu(1);
            x(type == 1 & typeold ==0) = x(type==1 & typeold ==0) + inj(type(typeold==0)==1) - kappau;

        end
        
        
          if t > Tinit
             
              typesim(:,t-Tinit) = type;
              agesim(:,t-Tinit)  = age;
              agemsim(:,t-Tinit) = agem;
              asim(:,t-Tinit)    = a;
              eesim(:,t-Tinit)   = e;
              xsim(:,t-Tinit)    = x + kappau*(type == 1 & typeold == 0) + kappap*(type == 2 & typeold == 1);    % convert back to savings prior to paying fixed cost
              xsim(type == 1 & typeold == 0,t-Tinit) = xsim(type == 1 & typeold == 0,t-Tinit) - inj(type(typeold==0)==1);
              typeoldsim(:,t-Tinit) = typeold;

          end
    end
end

a    = asim';
x    = xsim';
age  = agesim';
agem = agemsim';
e    = eesim';
type = typeoldsim';     % beginning of period type   
ru   = reshape(funeval(cru,fspaceu,[a(:), e(:)]), size(a,1), size(a,2));
rp   = reshape(funeval(crp,fspacep,[a(:), e(:)]), size(a,1), size(a,2));
rr   = ru.*(type==1) + rp.*(type==2) + r.*(type==0);
e    = egrid(e);


Lp   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e + phip);
Kp   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(e + phip);
Yp   = exp(e+phip).^(1-eta).*(Lp.^alpha.*Kp.^(1-alpha)).^eta;
Dp   = Yp - W.*Lp - (r+delta).*Kp;  

Lu   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e + phiu);
Ku   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(e + phiu);
Yu   = exp(e+phiu).^(1-eta).*(Lu.^alpha.*Ku.^(1-alpha)).^eta;
Du   = Yu - W.*Lu - (r+delta).*Ku;  

Lt   = (eta/W)^(1/(1-eta))*exp(e);
Yt   = (eta/W)^(eta/(1-eta))*exp(e);
Dt   = Yt - W.*Lt;

Y    = Yp.*(type==2) + Yu.*(type==1) + Yt.*(type==0);
K    = Kp.*(type==2) + Ku.*(type==1);
L    = Lp.*(type==2) + Lu.*(type==1) + Lt.*(type==0);
D    = Dp.*(type==2) + Du.*(type==1) + Dt.*(type==0);

cons = (1-theta*xai)*D + (1+r)*a - x;                    % consumption, (recall x is total savings)

pu   = reshape(funeval(cu(:,3), fspaceu, [a(:), e(:)]), size(a,1), size(a,2));
pp   = reshape(funeval(cp(:,3), fspacep, [a(:), e(:)]), size(a,1), size(a,2));

Eq     = 1/(1+r)*theta*xai*(D + pu.*(type==1) + pp.*(type==2));



% Frictionless allocations

Lpf   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e + phip);
Kpf   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(e + phip);
Ypf   = exp(e+phip).^(1-eta).*(Lpf.^alpha.*Kpf.^(1-alpha)).^eta;

Luf   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e + phiu);
Kuf   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(e + phiu);
Yuf   = exp(e+phiu).^(1-eta).*(Luf.^alpha.*Kuf.^(1-alpha)).^eta;

Ltf   = (eta/W)^(1/(1-eta))*exp(e);
Ytf   = (eta/W)^(eta/(1-eta))*exp(e);

Yf    = Ypf.*(type==2) + Yuf.*(type==1) + Ytf.*(type==0);
Kf    = Kpf.*(type==2) + Kuf.*(type==1);
Lf    = Lpf.*(type==2) + Luf.*(type==1) + Ltf.*(type==0);

A = a;

Debt   = K - A;
Equity = K + kappau.*(type==1) + (kappap+kappau).*(type==2) - Debt; 

% three ways to compute TFP losses in Modern Sector (TFP calculation only
% up to a scalar of normalization, may differ from number in ergodic.m)

et = e + phiu.*(type==1) + phip.*(type==2);   % total efficiency (exogenous and endogenous)

TFP      = mean(Y(type>0))./(mean(L(type>0)).^alpha*mean(K(type>0)).^(1-alpha))^eta;
TFPbest  = mean(Yf(type>0))./(mean(Lf(type>0)).^alpha*mean(Kf(type>0)).^(1-alpha))^eta;
TFP2best = mean(exp(et(type>0))).^(1-eta);

wedgel = (rr+delta).^(-(1-alpha)*eta/(1-eta)).*mean(exp(et(type>0)))./mean((rr(type>0)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(et(type>0)));
wedgek = (rr+delta).^((alpha*eta-1)/(1-eta)).*mean(exp(et(type>0)))./mean((rr(type>0)+delta).^((alpha*eta-1)/(1-eta)).*exp(et(type>0)));

TFPratio = mean(wedgel(type>0).^(alpha*eta).*wedgek(type>0).^((1-alpha)*eta).*exp(et(type>0)))./mean(exp(et(type>0)));

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
Kc = vec(K(2:end,:)) + kappau*0;
Kp = vec(K(1:end-1,:)) + kappau*0;

flag = vec(agem(2:end,:)>1); 

dY  = vec(log(Yc(flag)./Yp(flag))); 
dL  = vec(log(Lc(flag)./Lp(flag)));
dK  = vec(log(Kc(flag)./Kp(flag)));

% firm-level 'finance' statistics

YK3 = mean(log(Y(agem>10&type>0)./(K(agem>10&type>0)+kappau*0)));
YK1 = mean(log(Y(agem<=5&type>0)./(K(agem<=5&type>0)+kappau*0)))-YK3;
YK2 = mean(log(Y(agem>5&agem<=10&type>0)./(K(agem>5&agem<=10&type>0)+kappau*0)))-YK3;

YKa3 = log(sum(Y(agem>10&type>0))./sum(kappau*0 + K(agem>10&type>0)));
YKa1 = log(sum(Y(agem<=5&type>0))./sum(kappau*0 + K(agem<=5&type>0))) - YKa3;
YKa2 = log(sum(Y(agem>5&agem<=10&type>0))./sum(kappau*0 + K(agem>5&agem<=10&type>0))) - YKa3;

LK3 = mean(log(L(agem>10&type>0)./(K(agem>10&type>0)+kappau*0)));
LK1 = mean(log(L(agem<=5&type>0)./(K(agem<=5&type>0)+kappau*0)))-LK3;
LK2 = mean(log(L(agem>5&agem<=10&type>0)./(K(agem>5&agem<=10&type>0)+kappau*0)))-LK3;

LKa3 = log(sum(L(agem>10&type>0))./sum(kappau*0 + K(agem>10&type>0)));
LKa1 = log(sum(L(agem<=5&type>0))./sum(kappau*0 + K(agem<=5&type>0))) - LKa3;
LKa2 = log(sum(L(agem>5&agem<=10&type>0))./sum(kappau*0 + K(agem>5&agem<=10&type>0))) - LKa3;

Yc    = vec(Y(2:end,:));
Yp    = vec(Y(1:end-1,:));
Lc    = vec(L(2:end,:));
Lp    = vec(L(1:end-1,:));
Kc    = vec(K(2:end,:)) + kappau*0;
Kp    = vec(K(1:end-1,:)) + kappau*0;

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

if findvarz

varo  = var(et(type>0));
varz  = (0.0: .01: 2)';
varY  = var(log(Y(type>0)));
sY   = sqrt(varY + varz);
ac1   = (acorY.*varY  + varz)./(varY + varz);
ac3   = (acorY3.*varY + varz)./(varY + varz);
ac5   = (acorY5.*varY + varz)./(varY + varz);

obj = ((sY/sYdata - 1).^2 + (ac1/ac1Ydata -1).^2 + (ac3/ac3Ydata - 1).^2 + (ac5/ac5Ydata-1).^2);
[~, j] = min(obj); 

varz = varz(j);
sY  = sY(j); 
ac1  = ac1(j);
ac3  = ac3(j);
ac5  = ac5(j);

 save varz varz

else
    load varz
    
varo  = var(et(type>0));
varY  = var(log(Y(type>0)));
sdY   = sqrt(varY + varz);
ac1   = (acorY.*varY  + varz)./(varY + varz);
ac3   = (acorY3.*varY + varz)./(varY + varz);
ac5   = (acorY5.*varY + varz)./(varY + varz);
sY    = sqrt(varY + varz);
end
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

disp('Moments');

fprintf('\n');
disp('                   1st col: model,      2nd col: data');
fprintf('\n');
fprintf('s.d.  dlog(Y)                 = %9.2f  %9.2f \n',       [std(dY(:)),           sdYdata]);
fprintf('s.d.   log(Y)                 = %9.2f  %9.2f \n',       [sY,                   sYdata]);
fprintf('\n');
fprintf('cor y y1                      = %9.2f  %9.2f \n',       [ac1,                  ac1Ydata]);
fprintf('cor y y3                      = %9.2f  %9.2f \n',       [ac3,                  ac3Ydata]);
fprintf('cor y y5                      = %9.2f  %9.2f \n',       [ac5,                  ac5Ydata]);

fprintf('\n');
fprintf('aggregate Debt to GDP         = %9.2f  %9.2f \n',       [sum(Debt(type>0))./sum(Y(type>0)),       DtoY]); 
fprintf('average Debt to GDP |Debt>0   = %9.2f  %9.2f \n',       [sum(Debt(Debt>0 & type>0))./sum(Y(Debt>0&type>0)),       DtoY]); 

fprintf('aggregate posit. Debt to GDP  = %9.2f  %9.2f \n',       [sum(Debt(Debt>0 & type>0))./sum(Y(type>0)),       DtoY]); 

fprintf('Market capitaliz to GDP       = %9.2f  %9.2f \n',       [sum(Eq(type>0))/sum(Y(type>0)),          EtoY]);   
fprintf('\n');
finvest = mean(type(:)==1)*kappau*(mu-1) + mean(type(:)==2)*kappap*(mu-1);

fprintf('Fixed Inv  to total Y, perc   = %9.1f  %9.1f \n',       [finvest*numel(Y)/sum(Y(agem>=1))*100, 4.6]);
fprintf('\n');

fprintf('Variance exog component       = %9.2f \n',       varz);
% 
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf(' JIWOON, IGNORE THE REST OF THESE MOMENTS\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('\n');
% 
% 
% disp('Equilibrium Conditions')
% 
% fprintf('\n');
% 
% fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(K(:)),   mean(A(:)) + Aw - mean(Eq(:))]);
% fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(L(:)),   Lbar]);
% fprintf('\n');
% 
% fprintf('\n');
% disp('Model Implications')
% fprintf('\n');
% fprintf('Fraction productive           = %9.3f \n',          mean(type(:)==2));
% fprintf('Fraction unproductive         = %9.3f \n',          mean(type(:)==1));
% 
% fprintf('Fraction output in modern     = %9.3f \n',          sum(Y(type>0))/sum(Y(:)));
% fprintf('Fraction labor in modern      = %9.3f \n',          sum(L(type>0))/sum(L(:)));
% fprintf('\n');
% fprintf('Fract. var. perm comp.        = %9.2f \n',          varz/(varz+varo));
% 
% fprintf('\n');
% 
% Da = mean(K(:)) - mean(A(:)) - (Aw - mean(Eq(:))) ;       % Aggregate debt position
% 
% finvest = mean(type(:)==1)*kappau*(mu-1) + mean(type(:)==2)*kappap*(mu-1);
% 
% fprintf('Output                        = %9.3f \n',          mean(Y(:))/mean(L(:)));
% fprintf('Consumption                   = %9.3f \n',          (mean(Y(:)) - (delta+mu-1)*mean(K(:)) - (1+r-mu)*Da  - finvest)/mean(L(:)));
% fprintf('Investment                    = %9.3f \n',          ((delta + mu - 1)*mean(K(:)) + finvest)/mean(L(:)));
% fprintf('Fraction Fixed Investment     = %9.3f \n',          finvest/ ((delta + mu - 1)*mean(K(:)) + finvest));
% fprintf('TFP Traditional               = %9.3f \n',          sum(Y(type==0))/sum(L(type==0)).^(eta)*numel(Y)^(eta-1));   % express in per capita units, comparable to ergodic
% fprintf('TFP Modern                    = %9.3f \n',          sum(Y(type>0))./(sum(L(type>0)).^alpha*sum(K(type>0)).^(1-alpha))^eta*numel(Y)^(eta-1)); % express in per capita units, comparable to ergodic
% fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPbest/TFP)*100);


fprintf('\n');


ract = rr(type>0);

fprintf('\n');
fprintf('Fraction constrained          = %10.2f \n', mean(   ract > r+.00001));
fprintf('Median mu if constrained      = %10.2f \n', median( ract(ract > r+.00001)-r));
fprintf('Iqr mu if constrained         = %10.2f \n', iqr(    ract(ract > r+.00001)-r));
fprintf('90 pctile mu if constrained   = %10.2f \n', prctile(ract(ract > r+.00001)-r, 90));
fprintf('99 pctile mu if constrained   = %10.2f \n', prctile(ract(ract > r+.00001)-r, 99));
fprintf('Aggregate mu                  = %10.3f \n', sum(rr(type>0).*Y(type>0))./sum(Y(type>0))-r);


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
fprintf('fraction labor 1-5            = %9.2f  %9.2f \n',       [sum(L(agem(:)>=1 & agem(:)< 6))/sum(L(agem(:)>=1)),             0.20]);
fprintf('fraction labor 6-10           = %9.2f  %9.2f \n',       [sum(L(agem(:)>=6 & agem(:)< 11))/sum(L(agem(:)>=1)),            0.20]);
fprintf('fraction capital 1-5          = %9.2f  %9.2f \n',       [sum(K(agem(:)>=1 & agem(:)< 6))/sum(K(agem(:)>=1)),             0.20]);
fprintf('fraction capital 6-10         = %9.2f  %9.2f \n',       [sum(K(agem(:)>=6 & agem(:)< 11))/sum(K(agem(:)>=1)),            0.20]);

fprintf('\n');


fprintf('aggr DlogY 1-5                = %9.2f  %9.2f \n',     [dYa1-dYa3,     0.199]);
fprintf('aggr DlogY 6-10               = %9.2f  %9.2f \n',     [dYa2-dYa3,     0.089]);
fprintf('aggr DlogL 1-5                = %9.2f  %9.2f \n',     [dLa1-dLa3,     0.162]);
fprintf('aggr DlogL 6-10               = %9.2f  %9.2f \n',     [dLa2-dLa3,     0.068]);
fprintf('aggr DlogK 1-5                = %9.2f  %9.2f \n',     [dKa1-dKa3,     0.249]);
fprintf('aggr DlogK 6-10               = %9.2f  %9.2f \n',     [dKa2-dKa3,     0.040]);
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
sdL   = sqrt(varL + varz);
acL1  = (acorL.*varL  + varz)./(varL + varz);
acL3  = (acorL3.*varL + varz)./(varL + varz);
acL5  = (acorL5.*varL + varz)./(varL + varz);

varK  = var(log(kappau*0 + K(type>0)));
sdK   = sqrt(varK + varz);
acK1  = (acorK.*varK  + varz)./(varK + varz);
acK3  = (acorK3.*varK + varz)./(varK + varz);
acK5  = (acorK5.*varK + varz)./(varK + varz);

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

ffixed = kappau/(kappau + mean(K(agem==1)));

fprintf('Fixed to total K age = 1      = %9.2f  %9.2f \n',         [ffixed, 0.46]);

% Some implications that can test with data

fprintf('\n');
fprintf('\n');

fprintf('s.d.  dlog(Y)                 = %9.2f \n',     std(dY(:)));
fprintf('s.d.  dlog(K)                 = %9.2f \n',     std(dK(:)));
fprintf('\n');

flag = vec(agem(2:end,:)>1);  % everyone

Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));

Kc = vec(K(2:end,:));
Kp = vec(K(1:end-1,:));

ec = vec(e(2:end,:));
ep = vec(e(1:end-1,:));

dY  = log(Yc(flag)./Yp(flag)); 
dK  = log(Kc(flag)./Kp(flag));
de  = ec(flag) - ep(flag);


yy = dK; xx = [dY, ones(size(yy))];
bb = inv(xx'*xx)*xx'*yy; 

fprintf('elast dK to dY                = %9.2f \n',     bb(1));

yy = dY; xx = [de, ones(size(yy))];
bb = inv(xx'*xx)*xx'*yy; 

fprintf('elast dY to de                = %9.2f \n',     bb(1));

yy = dK; xx = [de, ones(size(yy))];
bb = inv(xx'*xx)*xx'*yy; 

fprintf('elast dK to de                = %9.2f \n',     bb(1));
fprintf('\n');

fprintf('DlogY 1-5 vs. 11 +            = %9.2f \n',     dYa1-dYa3);
fprintf('DlogK 1-5                     = %9.2f \n',     dKa1-dKa3);
fprintf('\n');
fprintf('mean  YK 1-5                  = %9.2f \n',     YK1);
fprintf('\n');

YK = log(Y(agem>=1)./K(agem>=1)); 
ee  = e(agem>=1); 

fprintf('var  YK                       = %9.3f \n',    var(YK));

% run regression of YK on age dummies

age = agem(agem>=1); 

yy = YK; xx = dummyvar(age);
bb = inv(xx'*xx)*xx'*yy; 

YKage = xx*bb; 
YKres = YK - YKage;  % not explained by age

fprintf('var  YK due to age            = %9.3f \n',    var(YKage));
fprintf('var  YK resid. 1-5            = %9.3f \n',    var(YKres(age<=5)));
fprintf('var  YK resid. 11+            = %9.3f \n',    var(YKres(age>10)));

TFPloss = mean(exp(ee - (1-alpha)*eta/(1-eta)*YK))^(1-alpha*eta)/...
          mean(exp(ee + (alpha*eta-1)/(1-eta)*YK))^((1-alpha)*eta)/mean(exp(ee))^(1-eta);

TFPlossage = mean(exp(ee - (1-alpha)*eta/(1-eta)*YKage))^(1-alpha*eta)/...
          mean(exp(ee + (alpha*eta-1)/(1-eta)*YKage))^((1-alpha)*eta)/mean(exp(ee))^(1-eta);

TFPlossy  = mean(exp(ee(age<=5) - (1-alpha)*eta/(1-eta)*YKres(age<=5)))^(1-alpha*eta)/...
          mean(exp(ee(age<=5) + (alpha*eta-1)/(1-eta)*YKres(age<=5)))^((1-alpha)*eta)/mean(exp(ee(age<=5)))^(1-eta);
 
TFPlosso  = mean(exp(ee(age>10) - (1-alpha)*eta/(1-eta)*YKres(age>10)))^(1-alpha*eta)/...
          mean(exp(ee(age>10) + (alpha*eta-1)/(1-eta)*YKres(age>10)))^((1-alpha)*eta)/mean(exp(ee(age>10)))^(1-eta);
      
      
fprintf('\n');
fprintf('\n');


fprintf('TFP loss                      = %9.3f \n',    -log(TFPloss)*100);
fprintf('TFP loss due to age           = %9.3f \n',    -log(TFPlossage)*100);
fprintf('TFP loss 1-5                  = %9.3f \n',    -log(TFPlossy)*100);
fprintf('TFP loss 11+                  = %9.3f \n',    -log(TFPlosso)*100);

fprintf('\n');
fprintf('\n');

flag = vec(agem(2:end,:)>1);  % everyone

Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));

Kc = vec(K(2:end,:));
Kp = vec(K(1:end-1,:));

ec = vec(e(2:end,:));
ep = vec(e(1:end-1,:));

dY  = log(Yc(flag)./Yp(flag)); 
dK  = log(Kc(flag)./Kp(flag));
de  = ec(flag) - ep(flag);

YK = log(Y(2:end,:)./K(2:end,:)); 
YK = YK(flag); 

fprintf('mean Y/K top 10p - low 10p dY        = %9.3f \n',    mean(YK(dY>prctile(dY,90))) - mean(YK(dY<prctile(dY,10))));
fprintf('mean Y/K top 10p - low 10p dK        = %9.3f \n',    mean(YK(dK>prctile(dK,90))) - mean(YK(dK<prctile(dK,10))));
fprintf('mean Y/K top 10p - low 10p de        = %9.3f \n',    mean(YK(de>prctile(de,90))) - mean(YK(de<prctile(de,10))));
fprintf('\n');
fprintf('var Y/K top 10p dY        = %9.3f \n',    var(YK(dY>prctile(dY,90))));
fprintf('var Y/K botton 10p dY     = %9.3f \n',    var(YK(dY<prctile(dY,10))));

fprintf('var Y/K top 10p dK        = %9.3f \n',    var(YK(dK>prctile(dK,90))));
fprintf('var Y/K botton 10p dK     = %9.3f \n',    var(YK(dK<prctile(dK,10))));

fprintf('var Y/K top 10p de        = %9.3f \n',    var(YK(de>prctile(de,90))));
fprintf('var Y/K botton 10p de     = %9.3f \n',    var(YK(de<prctile(de,10))));
fprintf('\n');
fprintf('\n');


flag = vec(agem(2:end,:)>1 & agem(2:end,:)<=5);  % young
dY  = log(Yc(flag)./Yp(flag)); 

fprintf('dY young: s.d., 25th, 50th, 75th        = %9.3f %9.3f %9.3f %9.3f \n',   [std(dY), prctile(dY,25), prctile(dY,50), prctile(dY, 75)]);


flag = vec(agem(2:end,:)>10);  % young
dY  = log(Yc(flag)./Yp(flag)); 

fprintf('dY old: s.d., 25th, 50th, 75th          = %9.3f %9.3f %9.3f %9.3f \n',   [std(dY), prctile(dY,25), prctile(dY,50), prctile(dY, 75)]);
