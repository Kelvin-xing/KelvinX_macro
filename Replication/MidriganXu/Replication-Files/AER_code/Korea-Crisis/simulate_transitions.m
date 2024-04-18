Tinit = 100;       % periods to initialize distribution of assets
Nf    = 500000;    % # producers

randn('state', 100);
rand('state',  100);

eesim    = zeros(Nf,T);       % transitory component
typesim  = zeros(Nf,T);       % end of period type: 2 = productive, 1 = unproductive, 0 = traditional 
typeoldsim = zeros(Nf,T);     % beginning of period type
agesim   = zeros(Nf,T);       % age
agemsim  = zeros(Nf,T);       % age in modern sector
asim     = zeros(Nf,T);       % beginning of period assets
xsim     = zeros(Nf,T);        
a        = ones(Nf,1);

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

for t = 1 : Tinit + T

    if t==1
          
        e      = esim(:,t);

       % run this to make sure same shocks
                
    else

        jump   = jump(randperm(Nf));
       
        e      = esim(:,t);                  % "exit" is random, so ok to use incumbent's e as draw from ergodic

        % Tinit is date for which we save the state
        % Tinit + 1 is using pre-crisis coefficients
        
        if t > Tinit
           
            if t == Tinit + 1
                
         load c_pre; 
                
         load sim_steady
         type = sim_steady(:,1);
         age  = sim_steady(:,2);
         agem = sim_steady(:,3);
         x    = sim_steady(:,4);
       
         
        type   = type.*(~jump);

        age    = age.*(~jump) + 1;
        agem   = (agem + 1).*(type>=1);
        
        x      = x.*(~jump) + smintold(1)*(jump);
 
  typeold    = type;   % save to determine who pays the fixed cost
        
        if any(typeold == 1)
            
            state      = [x(typeold == 1), e(typeold == 1)];
            a(typeold == 1) = state(:,1);
            x(typeold == 1) = max(funeval(cxu, fspaceuold, state), sminuold(1)); 

            type(typeold == 1) = type(typeold == 1);

        end
        
        if any(typeold == 0)
                        
            state         = [x(typeold == 0), e(typeold==0)];
            a(typeold==0) = state(:,1);
            x(typeold==0) = max(funeval(cxt, fspacetold, state), smintold(1)); 
           
            Phi           = funbas(fspacetold, [x(typeold==0), e(typeold==0)]);
            
            type(typeold==0) = Phi*cst(:,1)>Phi*cst(:,2) & x(typeold==0) - kappau > sminuold(1);
            x(type == 1 & typeold ==0) = x(type==1 & typeold ==0) - kappau;

        end
         
            else
                            
    cxu = cxu_irf(:,:,t - Tinit - 1); 
    cxt = cxt_irf(:,:,t - Tinit - 1);
    cst = cst_irf(:,:,t - Tinit - 1);
        
        type   = type.*(~jump);

        age    = age.*(~jump) + 1;
        agem   = (agem + 1).*(type>=1);
        
        x      = x.*(~jump) + smint(1)*(jump);
 
  typeold    = type;   % save to determine who pays the fixed cost
        
        if any(typeold == 1)
            
            state      = [x(typeold == 1), e(typeold == 1)];
            a(typeold == 1) = state(:,1);
            x(typeold == 1) = max(funeval(cxu, fspaceu, state), sminu(1)); 

            type(typeold == 1) = type(typeold == 1);

        end
        
        if any(typeold == 0)
                        
            state         = [x(typeold == 0), e(typeold==0)];
            a(typeold==0) = state(:,1);
            x(typeold==0) = max(funeval(cxt, fspacet, state), smint(1)); 
           
            Phi           = funbas(fspacet, [x(typeold==0), e(typeold==0)]);
            
            type(typeold==0) = Phi*cst(:,1)>Phi*cst(:,2) & x(typeold==0) - kappau > sminu(1);
            x(type == 1 & typeold ==0) = x(type==1 & typeold ==0) - kappau;

        end
                
            end
              typesim(:,t-Tinit) = type;
              agesim(:,t-Tinit)  = age;
              agemsim(:,t-Tinit) = agem;
              asim(:,t-Tinit)    = a;
              eesim(:,t-Tinit)   = e;
              xsim(:,t-Tinit)    = x + kappau*(type == 1 & typeold == 0);    % convert back to savings prior to paying fixed cost
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

rr   = zeros(T,Nf);

for t = 1 : T
    
    if t == 1
        
        load c_pre   % loads cru from prior to crisis

ru   = funeval(cru,fspaceuold, [a(t,:)', e(t,:)'])';
rr(t,:)   = ru.*(type(t,:)==1) + r.*(type(t,:)==0);

    else
        
cru = cru_irf(:,:,t-1);

ru  = funeval(cru,fspaceu, [a(t,:)', e(t,:)'])';
rr(t,:)   = ru.*(type(t,:)==1) + r.*(type(t,:)==0);

    end    

end

e    = egrid(e);

Lu = zeros(T,Nf);
Ku = Lu;
Yu = Lu;
Du = Lu;
Lt = Lu;
Yt = Lu;
Dt = Lu;

Luf = Lu;
Kuf = Ku; 
Yuf = Yu;
Ltf = Lu; 
Ytf = Yu; 

for t = 1 : T

if t == 1
    
W      = 1.042530270930890;

else
    
W   = Wsave(t-1);

end
    
Lu(t,:)   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr(t,:)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e(t,:) + phiu);
Ku(t,:)   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr(t,:)+delta).^((alpha*eta-1)/(1-eta)).*exp(e(t,:) + phiu);
Yu(t,:)   = exp(e(t,:)+phiu).^(1-eta).*(Lu(t,:).^alpha.*Ku(t,:).^(1-alpha)).^eta;
Du(t,:)   = Yu(t,:) - W.*Lu(t,:) - (r+delta).*Ku(t,:);  

Lt(t,:)   = (eta/W)^(1/(1-eta))*exp(e(t,:));
Yt(t,:)   = (eta/W)^(eta/(1-eta))*exp(e(t,:));
Dt(t,:)   = Yt(t,:) - W.*Lt(t,:);

% Frictionless allocations

Luf(t,:)   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e(t,:) + phiu);
Kuf(t,:)   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(e(t,:) + phiu);
Yuf(t,:)   = exp(e(t,:)+phiu).^(1-eta).*(Luf(t,:).^alpha.*Kuf(t,:).^(1-alpha)).^eta;

Ltf(t,:)   = (eta/W)^(1/(1-eta))*exp(e(t,:));
Ytf(t,:)   = (eta/W)^(eta/(1-eta))*exp(e(t,:));

end

% Y    = Yu.*(type==1) + Yt.*(type==0);
% K    = Ku.*(type==1);
% L    = Lu.*(type==1) + Lt.*(type==0);
% D    = Du.*(type==1) + Dt.*(type==0);
% 
% cons = D + (1+r)*a - x;                    % consumption, (recall x is total savings)
% 
% Yf    = Yuf.*(type==1) + Ytf.*(type==0);
% Kf    = Kuf.*(type==1);
% Lf    = Luf.*(type==1) + Ltf.*(type==0);

Yua       = mean(Yu.*(type==1),2);
Kua       = mean(Ku.*(type==1),2);
Lua       = mean(Lu.*(type==1),2); 
Yta       = mean(Yt.*(type==0),2);
Lta       = mean(Lt.*(type==0),2);
Aua       = mean(a.*(type==1),2);
Ata       = mean(a.*(type==0),2);
Yufa      = mean(Yuf.*(type==1),2);
Kufa      = mean(Kuf.*(type==1),2);
Lufa      = mean(Luf.*(type==1),2);
TFPu      = Yua./(Lua.^alpha.*Kua.^(1-alpha)).^(eta);
TFPt      = Yta./Lta.^eta;

et = e + phiu;   % total efficiency (exogenous and endogenous)

TFPubest  = mean(exp(et).*(type==1),2).^(1-eta);

% TFP losses due to age in year of the crisis: 

% use average for each age

rdelta = zeros(max(agem(2,:)'), 1);
rrnew = zeros(Nf,1);

for i = 1 : max(agem(2,:))
    
    rdelta(i) = mean(rr(2,agem(2,:)==i & type(2,:)>0)); 
    rrnew(agem(2,:)==i) = rdelta(i);
    
end

et    = et(2,rrnew>0)';
rrnew = rrnew(rrnew>0); 

TFPlossage = mean(exp(et).*(rrnew+delta).^(-(1-alpha)*eta/(1-eta))).^(1-alpha*eta)/...
    mean(exp(et).*(rrnew+delta).^((alpha*eta-1)/(1-eta))).^((1-alpha)*eta)/mean(exp(et)).^(1-eta);

fprintf('TFP due to age                = %9.2f \n',     -log(TFPlossage)*100);


num       = mean(type==1,2);       % measure modern sector

% compute statistics year by year. Then add statistics prior to crisis

stdY    = zeros(T,1);
stdK    = zeros(T,1);
stdL    = zeros(T,1);
stdYK   = zeros(T,1);
stdLK   = zeros(T,1);
 
stdYKyoung   = zeros(T,1);
stdYKold   = zeros(T,1);



YKyoung = zeros(T,1);
LKyoung = zeros(T,1);
YKold   = zeros(T,1);
LKold   = zeros(T,1);

dYyoung = zeros(T,1);
dKyoung = zeros(T,1);
dLyoung = zeros(T,1);
dYold   = zeros(T,1);
dKold   = zeros(T,1);
dLold   = zeros(T,1);

eYK      = zeros(T,1);    % sensitivity of YK to delta(Y)
eYKyoung = zeros(T,1);
eYKold   = zeros(T,1);

eYKe      = zeros(T,1);   % sensitivity of YK to delta(e)
eYKeyoung = zeros(T,1);
eYKeold   = zeros(T,1);

eIK       = zeros(T,1);    % sensitivity of IK to YK


for t = 1:T
   
    YK = log(Yu(t,agem(t,:)>=1)./Ku(t,agem(t,:)>=1));
    LK = log(Lu(t,agem(t,:)>=1)./Ku(t,agem(t,:)>=1));
    
    stdYK(t) = std(YK);
    stdLK(t) = std(LK);
    
    
    YK = log(Yu(t,agem(t,:)>=1&agem(t,:)<=5)./Ku(t,agem(t,:)>=1&agem(t,:)<=5));
    
    stdYKyoung(t) = std(YK); 
    
    YK = log(Yu(t,agem(t,:)>=11)./Ku(t,agem(t,:)>=11));

    stdYKold(t) = std(YK); 

    
    YKyoung(t) = log(sum(Yu(t,agem(t,:)>=1&agem(t,:)<=5))./sum(Ku(t,agem(t,:)>=1&agem(t,:)<=5)));
    LKyoung(t) = log(sum(Lu(t,agem(t,:)>=1&agem(t,:)<=5))./sum(Ku(t,agem(t,:)>=1&agem(t,:)<=5)));
    YKold(t) = log(sum(Yu(t,agem(t,:)>10))./sum(Ku(t,agem(t,:)>=10)));
    LKold(t) = log(sum(Lu(t,agem(t,:)>10))./sum(Ku(t,agem(t,:)>=10)));
    
    if t > 1
    
    dY = log(Yu(t,agem(t,:)>1)./Yu(t-1,agem(t,:)>1));
    dK = log(Ku(t,agem(t,:)>1)./Ku(t-1,agem(t,:)>1));
    dL = log(Lu(t,agem(t,:)>1)./Lu(t-1,agem(t,:)>1));
    
    stdY(t) = std(dY);
    stdK(t) = std(dK);
    stdL(t) = std(dL);
    
    dYyoung(t) = log(sum(Yu(t,agem(t,:)>=2&agem(t,:)<=5))./sum(Yu(t-1,agem(t,:)>=2&agem(t,:)<=5)));
    dYold(t)   = log(sum(Yu(t,agem(t,:)>10))./sum(Yu(t-1,agem(t,:)>10)));
    dLyoung(t) = log(sum(Lu(t,agem(t,:)>=2&agem(t,:)<=5))./sum(Lu(t-1,agem(t,:)>=2&agem(t,:)<=5)));
    dLold(t)   = log(sum(Lu(t,agem(t,:)>10))./sum(Lu(t-1,agem(t,:)>10)));
    dKyoung(t) = log(sum(Ku(t,agem(t,:)>=2&agem(t,:)<=5))./sum(Ku(t-1,agem(t,:)>=2&agem(t,:)<=5)));
    dKold(t)   = log(sum(Ku(t,agem(t,:)>10))./sum(Ku(t-1,agem(t,:)>10)));
 
    % elasticity Y/K to delta Y
    
    cyk = log(Yu(t,agem(t,:)>1)./Ku(t,agem(t,:)>1))';
    lyk = log(Yu(t-1,agem(t,:)>1)./Ku(t-1,agem(t,:)>1))';
    dy  = log(Yu(t,agem(t,:)>1)./Yu(t-1,agem(t,:)>1))';

    xx = [ones(size(cyk)), dy, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYK(t) = beta(2);
      
    % elasticity Y/K to delta Y for young
    
    cyk = log(Yu(t,agem(t,:)>1&agem(t,:)<=5)./Ku(t,agem(t,:)>1&agem(t,:)<=5))';
    lyk = log(Yu(t-1,agem(t,:)>1&agem(t,:)<=5)./Ku(t-1,agem(t,:)>1&agem(t,:)<=5))';
    dy  = log(Yu(t,agem(t,:)>1&agem(t,:)<=5)./Yu(t-1,agem(t,:)>1&agem(t,:)<=5))';

    xx = [ones(size(cyk)), dy, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYKyoung(t) = beta(2);
    
    % elasticity Y/K to delta Y for young

    cyk = log(Yu(t,agem(t,:)>10)./Ku(t,agem(t,:)>10))';
    lyk = log(Yu(t-1,agem(t,:)>10)./Ku(t-1,agem(t,:)>10))';
    dy  = log(Yu(t,agem(t,:)>10)./Yu(t-1,agem(t,:)>10))';

    xx = [ones(size(cyk)), dy, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYKold(t) = beta(2);
    
     % elasticity Y/K to delta e
    
    cyk = log(Yu(t,agem(t,:)>1)./Ku(t,agem(t,:)>1))';
    lyk = log(Yu(t-1,agem(t,:)>1)./Ku(t-1,agem(t,:)>1))';
    de  = e(t,agem(t,:)>1)' - e(t-1,agem(t,:)>1)';

    xx = [ones(size(cyk)), de, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYKe(t) = beta(2);
      
    % elasticity Y/K to delta e for young
    
    cyk = log(Yu(t,agem(t,:)>1&agem(t,:)<=5)./Ku(t,agem(t,:)>1&agem(t,:)<=5))';
    lyk = log(Yu(t-1,agem(t,:)>1&agem(t,:)<=5)./Ku(t-1,agem(t,:)>1&agem(t,:)<=5))';
    de  = e(t,agem(t,:)>1&agem(t,:)<=5)' - e(t-1,agem(t,:)>1&agem(t,:)<=5)';

    xx = [ones(size(cyk)), de, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYKeyoung(t) = beta(2);
    
    % elasticity Y/K to delta e for young
    
    cyk = log(Yu(t,agem(t,:)>10)./Ku(t,agem(t,:)>10))';
    lyk = log(Yu(t-1,agem(t,:)>10)./Ku(t-1,agem(t,:)>10))';
    de  = e(t,agem(t,:)>10)' - e(t-1,agem(t,:)>10)';

    xx = [ones(size(cyk)), de, lyk];
    yy = cyk;

    beta = inv(xx'*xx)*xx'*yy;
    
    eYKeold(t) = beta(2);
    
    end
    
    if t < T
        
       ik = (Ku(t+1,agem(t+1,:)>1) -  Ku(t,agem(t+1,:)>1))'./Ku(t,agem(t+1,:)>1)';
       yk = Yu(t,agem(t+1,:)>1)'./Ku(t,agem(t+1,:)>1)';
        
        xx = [ones(size(ik)), log(yk)];
        yy = ik;
        
        beta = inv(xx'*xx)*xx'*yy;

        
        eIK(t) = beta(2);
       
    end
    
end

close all
set(gcf,'DefaultLineLineWidth',2);
subplot(3,2,1)
hold on
plot((1:1:10)',[stdY(2:11), stdK(2:11)]);
legend('s.d. delta Y', 's.d. delta K')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)

subplot(3,2,2)
hold on
plot((0:1:10)', stdYK(1:11))
legend('s.d. Y/K')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)

subplot(3,2,3)
hold on
plot((0:1:10)', [YKyoung(1:11), YKold(1:11)])
legend('YK young', 'YK old')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)

subplot(3,2,4)
hold on
plot((1:1:10)', [dYyoung(2:11), dYold(2:11)])
legend('dY young', 'dY old')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)

subplot(3,2,5)
hold on
plot((1:1:10)', [eYKyoung(2:11), eYKold(2:11)])
legend('elasticity YK to Y young', 'elasticity YK to Y old')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)

subplot(3,2,6)
hold on
plot((1:1:10)', [eYKeyoung(2:11), eYKeold(2:11)])
legend('elasticity YK to e young', 'elasticity YK to e old')
xlabel('years')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',2)
