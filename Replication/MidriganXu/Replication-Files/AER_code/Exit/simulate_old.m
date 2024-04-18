Tinit = 100;       % periods to initialize distribution of assets
T     = 50;        % periods to simulate allocations
Nf    = 50000;     % # producers

randn('state', 100);
rand('state',  100);

eesim    = zeros(Nf,T);       % transitory component
typesim  = zeros(Nf,T);       % 1 = entrepreneur, 0 = worker 
agesim   = zeros(Nf,T);
agemsim  = zeros(Nf,T);
asim     = zeros(Nf,T);

Pcum      = P*triu(ones(size(P)));   % cumulative ergodic distributon

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

type  = zeros(Nf,1);    % start everyone as worker 

for t = 1 : Tinit+T

    if t==1
         
        e      = esim(:,t);
        age    = ones(Nf,1);
        agem   = zeros(Nf,1);
        
        state  = [smint(1)*ones(Nf,1), e]; 
        a      = state(:,1);
        x      = funeval(cxt,fspacet,state);

        Phi    = funbas(fspacet,[x,e]);        % Evaluate value of switching/staying next period give current x
        
        type   = Phi*csw(:,1)>Phi*csw(:,2) & x - kappa > sminm(1);
        x(type==1) = x(type==1) - kappa; 
                
    else
      
        jump   = jump(randperm(Nf));
       
        e      = esim(:,t);                  % "exit" is random, so ok to use incumbent's e as draw from ergodic
        type   = type.*(~jump);

        age    = age.*(~jump) + 1;
        agem   = (agem + 1).*(type==1);
        
        x      = x.*(~jump) + smint(1)*(jump);
        
        if any(type == 1)
            
           state      = [x(type == 1), e(type==1)];
           a(type==1) = state(:,1);
           x(type==1) = funeval(cxm, fspacem, state);
            
        end
        
        if any(type == 0)
            
            typeold    = type;   % save to determine who pays the fixed cost
            
            state      = [x(type == 0), e(type==0)];
            a(type==0) = state(:,1);
            x(type==0) = funeval(cxt, fspacet, state);
           
            Phi        = funbas(fspacet, [x(type==0), e(type==0)]);

            type(type==0) = Phi*csw(:,1)>Phi*csw(:,2) & x(type==0) - kappa > sminm(1);
            x(type == 1 & typeold ==0) = x(type==1 & typeold ==0) - kappa;

        end
        
          if t > Tinit
             
              typesim(:,t-Tinit) = type;
              agesim(:,t-Tinit)  = age;
              agemsim(:,t-Tinit) = agem;
              asim(:,t-Tinit)    = a;
              eesim(:,t-Tinit)   = e;
              
          end
    end
end

a    = asim';
age  = agesim';
agem = agemsim';
e    = eesim';
rr   = reshape(funeval(cr,fspacer,[a(:), e(:)]), size(a,1), size(a,2));
e    = egrid(e);
type = agem>0;

% z  = nodeunif(Nf, eps^(1/3.35), 1-eps^(1/3.35));  
% z  = log(z.^(-1/sz));                              % Pareto

z  = nodeunif(Nf, eps, 1-eps);
u1 = z(randperm(Nf));
u2 = z(randperm(Nf));
z  = sqrt(-2*log(u1)).*cos(2*pi*u2);                 % Gaussian with st.dev. 1
z(1:floor(pz*Nf)) = z(1:floor(pz*Nf))*sz;
z(floor(pz*Nf)+1:end) = z(floor(pz*Nf)+1:end)*lz*sz;


if 0   
    
   yy = [sort(exp(z)), (Nf:-1:1)'./Nf]; loglog(yy(:,1),yy(:,2),'bo');
   hold on
   plot(yy(:,1),yy(:,1).^(-mu),'r--');   
   
end

z   = z(randperm(Nf));
z   = repmat(z',T,1);

L   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e+z);
K   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(e+z);
Y   = exp(e+z).^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D   = Y - W.*L - (r+delta).*K;  
    
Y   = Y.*(agem>0);
K   = K.*(agem>0);
L   = L.*(agem>0);

Lf  = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e+z);
Kf  = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(e+z);
Yf  = exp(e+z).^(1-eta).*(Lf.^alpha.*Kf.^(1-alpha)).^eta;

Df  = Yf - W.*Lf - (r+delta).*Kf;  

Yf  = Yf.*(agem>0);
Kf  = Kf.*(agem>0);
Lf  = Lf.*(agem>0);

A  = a.*exp(z);
Debt = K - A;
Equity = K + exp(z).*kappa.*(type>0) - Debt; 
DebtRatio  = zeros(T,Nf);
DebtRatio(Y>0)  = Debt(Y>0)./Equity(Y>0); 

% three ways to compute TFP losses

TFP     = mean(Y(Y>0))./(mean(L(L>0)).^alpha*mean(K(K>0)).^(1-alpha))^eta;
TFPbest = mean(Yf(Y>0))./(mean(Lf(Lf>0)).^alpha*mean(Kf(Kf>0)).^(1-alpha))^eta;
TFP2best = mean(exp(z(Y>0)+e(Y>0))).^(1-eta);

wedgel = (rr+delta).^(-(1-alpha)*eta/(1-eta)).*mean(exp(z(Y>0)+e(Y>0)))./mean((rr(Y>0)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(z(Y>0)+e(Y>0)));
wedgek = (rr+delta).^((alpha*eta-1)/(1-eta)).*mean(exp(z(Y>0)+e(Y>0)))./mean((rr(Y>0)+delta).^((alpha*eta-1)/(1-eta)).*exp(z(Y>0)+e(Y>0)));

TFPratio = mean(wedgel(Y>0).^(alpha*eta).*wedgek(Y>0).^((1-alpha)*eta).*exp(e(Y>0)+z(Y>0)))./mean(exp(z(Y>0)+e(Y>0)));

% autocorrelation Y 

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

% concentration statistics:  

data = sort(vec(Y(Y>0)),1);

fRtopY20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopY10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopY5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopY1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(Y(Y>0&agem<=5)),1);

fRtopYy20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopYy10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopYy5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopYy1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(Y(Y>0&agem>5&agem<=10)),1);

fRtopYm20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopYm10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopYm5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopYm1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(Y(Y>0 & agem>10)),1);

fRtopYo20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopYo10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopYo5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopYo1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(L(L>0)),1);

fRtopL20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopL10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopL5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopL1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);


data = sort(vec(L(L>0&agem<=5)),1);

fRtopLy20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopLy10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopLy5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopLy1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(L(L>0&agem>5&agem<=10)),1);

fRtopLm20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopLm10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopLm5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopLm1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(L(L>0 & agem>10)),1);

fRtopLo20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopLo10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopLo5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopLo1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(K(K>0)),1);

fRtopK20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopK10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopK5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopK1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(K(K>0&agem<=5)),1);

fRtopKy20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopKy10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopKy5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopKy1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(K(K>0&agem>5&agem<=10)),1);

fRtopKm20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopKm10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopKm5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopKm1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

data = sort(vec(K(K>0 & agem>10)),1);

fRtopKo20 = sum(data(end-floor(numel(data)*0.20)+1:end))/sum(data);
fRtopKo10 = sum(data(end-floor(numel(data)*0.10)+1:end))/sum(data);
fRtopKo5  = sum(data(end-floor(numel(data)*0.05)+1:end))/sum(data);
fRtopKo1  = sum(data(end-floor(numel(data)*0.01)+1:end))/sum(data);

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

YK  = vec(log(Y(Y>0)./K(Y>0)));  
LK  = vec(log(L(Y>0)./K(Y>0)));  


YK3 = mean(log(Y(agem>10&Y>0)./K(agem>10&Y>0)));
YK1 = mean(log(Y(agem<=5&Y>0)./K(agem<=5&Y>0)))-YK3;
YK2 = mean(log(Y(agem>5&agem<=10&Y>0)./K(agem>5&agem<=10&Y>0)))-YK3;

YKa3 = log(sum(Y(agem>10&Y>0))./sum(K(agem>10&Y>0)));
YKa1 = log(sum(Y(agem<=5&Y>0))./sum(K(agem<=5&Y>0))) - YKa3;
YKa2 = log(sum(Y(agem>5&agem<=10&Y>0))./sum(K(agem>5&agem<=10&Y>0))) - YKa3;

LK3 = mean(log(L(agem>10&Y>0)./K(agem>10&Y>0)));
LK1 = mean(log(L(agem<=5&Y>0)./K(agem<=5&Y>0)))-LK3;
LK2 = mean(log(L(agem>5&agem<=10&Y>0)./K(agem>5&agem<=10&Y>0)))-LK3;

LKa3 = log(sum(L(agem>10&Y>0))./sum(K(agem>10&Y>0)));
LKa1 = log(sum(L(agem<=5&Y>0))./sum(K(agem<=5&Y>0))) - LKa3;
LKa2 = log(sum(L(agem>5&agem<=10&Y>0))./sum(K(agem>5&agem<=10&Y>0))) - LKa3;

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

disp('Moments');
fprintf('\n');
disp('1st col: model, 2nd col: data');
fprintf('\n');
fprintf('s.d.  dlog(Y)                 = %9.2f  %9.2f \n',       [std(dY(:)),           0.56]);
fprintf('s.d.   log(Y)                 = %9.2f  %9.2f \n',       [std(log(Y(Y>0))),     1.29]);
fprintf('\n');
fprintf('cor y y1                      = %9.2f  %9.2f \n',       [acorY,                0.91]);
fprintf('cor y y3                      = %9.2f  %9.2f \n',       [acorY3,               0.87]);
fprintf('cor y y5                      = %9.2f  %9.2f \n',       [acorY5,               0.85]);
fprintf('\n');
fprintf('Share of Y of top 0.01        = %9.2f  %9.2f \n',       [fRtopY1,              0.54]);
fprintf('Share of Y of top 0.05        = %9.2f  %9.2f \n',       [fRtopY5,              0.73]);
fprintf('Share of Y of top 0.1         = %9.2f  %9.2f \n',       [fRtopY10,             0.81]);
fprintf('Share of Y of top 0.2         = %9.2f  %9.2f \n',       [fRtopY20,             0.88]);
fprintf('\n');
fprintf('aggr DlogY 1-5                = %9.2f  %9.2f \n',       [dYa1-dYa3,    0.140  - 0.056]);
fprintf('aggr DlogY 6-10               = %9.2f  %9.2f \n',       [dYa2-dYa3,    0.078  - 0.056]);
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('aggregate debt to GDP         = %9.1f  %9.1f \n',       [sum(Debt(Debt>0))./sum(Y(Y>0)),  1.2]);    
fprintf('aggregate Debt to Equity      = %9.1f  %9.1f \n',       [sum(Debt(Y>0))./sum(Equity(Y>0)),  3.01]);   

fprintf('\n');
fprintf('\n');
disp('Implications for TFP')
fprintf('\n');
fprintf('TFP losses  misall           = %10.2f \n', -100*(log(TFPratio)));

ract = rr(Y>0);
fprintf('\n');
fprintf('Fraction constrained         = %10.2f \n', mean(   ract > r + 0.0025));
fprintf('Median r if constrained      = %10.2f \n', median( ract(ract > r + 0.0025)));
fprintf('Iqr r if constrained         = %10.2f \n', iqr(    ract(ract > r + 0.0025)));
fprintf('90 pctile r if constrained   = %10.2f \n', prctile(ract(ract > r + 0.0025),90));
fprintf('95 pctile r if constrained   = %10.2f \n', prctile(ract(ract > r + 0.0025),95));
fprintf('99 pctile r if constrained   = %10.2f \n', prctile(ract(ract > r + 0.0025),99));
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

disp('Other moments')
fprintf('\n');

fprintf('s.d.   log(Y) young           = %9.2f  %9.2f \n',     [std(log(Y(agem<=5&Y>0))),  1.08]);
fprintf('s.d.   log(Y) old             = %9.2f  %9.2f \n',     [std(log(Y(agem>10&Y>0))),  1.62]);
fprintf('s.d.   log(L) young           = %9.2f  %9.2f \n',     [std(log(L(agem<=5&Y>0))),  1.21]);
fprintf('s.d.   log(L) old             = %9.2f  %9.2f \n',     [std(log(L(agem>10&Y>0))),  1.68]);
fprintf('s.d.   log(K) young           = %9.2f  %9.2f \n',     [std(log(K(agem<=5&Y>0))),  1.21]);
fprintf('s.d.   log(K) old             = %9.2f  %9.2f \n',     [std(log(K(agem>10&Y>0))),  1.68]);
fprintf('\n');

fprintf('mean DlogY 1-5                = %9.2f  %9.2f \n',     [dY1-dY3,    0.084 - 0.03]);
fprintf('mean DlogY 6-10               = %9.2f  %9.2f \n',     [dY2-dY3,    0.054 - 0.03]);
fprintf('mean DlogL 1-5                = %9.2f  %9.2f \n',     [dL1-dL3,    0.066 - 0.022]);
fprintf('mean DlogL 6-10               = %9.2f  %9.2f \n',     [dL2-dL3,    0.041 - 0.022]);
fprintf('mean DlogK 1-5                = %9.2f  %9.2f \n',     [dK1-dK3,    0.047 - 0.009]);
fprintf('mean DlogK 6-10               = %9.2f  %9.2f \n',     [dK2-dK3,    0.018 - 0.009]);
fprintf('\n');
fprintf('\n');
fprintf('aggr DlogY 1-5                = %9.2f  %9.2f \n',     [dYa1-dYa3,    0.140 - 0.056]);
fprintf('aggr DlogY 6-10               = %9.2f  %9.2f \n',     [dYa2-dYa3,    0.078 - 0.056]);
fprintf('aggr DlogL 1-5                = %9.2f  %9.2f \n',     [dLa1-dLa3,    0.103 - 0.032]);
fprintf('aggr DlogL 6-10               = %9.2f  %9.2f \n',     [dLa2-dLa3,    0.053 - 0.032]);
fprintf('aggr DlogK 1-5                = %9.2f  %9.2f \n',     [dKa1-dKa3,    0.092 - 0.077]);
fprintf('aggr DlogK 6-10               = %9.2f  %9.2f \n',     [dKa2-dKa3,    0.040 - 0.077]);
fprintf('\n');
fprintf('mean  YK 1-5                  = %9.2f  %9.2f \n',     [YK1,     -0.03]);
fprintf('mean  YK 6-10                 = %9.2f  %9.2f \n',     [YK2,     0]);
fprintf('mean  LK 1-5                  = %9.2f  %9.2f \n',     [LK1,    -0.03]);
fprintf('mean  LK 6-10                 = %9.2f  %9.2f \n',     [LK2,    -0.01]);
fprintf('\n');
fprintf('aggr  YK 1-5                  = %9.2f  %9.2f \n',     [YKa1,    -0.16]);
fprintf('aggr  YK 6-10                 = %9.2f  %9.2f \n',     [YKa2,    -0.01]);
fprintf('aggr  LK 1-5                  = %9.2f  %9.2f \n',     [LKa1,    -0.11]);
fprintf('aggr  LK 6-10                 = %9.2f  %9.2f \n',     [LKa2,    -0.02]);
fprintf('\n');
fprintf('\n');
fprintf('s.d.  dlog(L)                 = %9.2f  %9.2f \n',     [std(dL),        0.44]);
fprintf('s.d.   log(L)                 = %9.2f  %9.2f \n',     [std(log(L(Y>0))),    1.19]);
fprintf('\n');
fprintf('s.d.  dlog(K)                 = %9.2f  %9.2f \n',     [std(dK),        0.40]);
fprintf('s.d.   log(K)                 = %9.2f  %9.2f \n',     [std(log(K(Y>0))),    1.39]);
fprintf('\n');
fprintf('cor L L1                      = %9.2f  %9.2f \n',     [acorL,             0.94]);
fprintf('cor L L3                      = %9.2f  %9.2f \n',     [acorL3,            0.90]);
fprintf('cor L L5                      = %9.2f  %9.2f \n',     [acorL5,            0.87]);
fprintf('\n');
fprintf('cor K K1                      = %9.2f  %9.2f \n',     [acorK,             0.96]);
fprintf('cor K K3                      = %9.2f  %9.2f \n',     [acorK3,            0.92]);
fprintf('cor K K5                      = %9.2f  %9.2f \n',     [acorK5,            0.89]);
fprintf('\n');
fprintf('Share of L of top 0.01        = %9.2f  %9.2f \n',     [fRtopL1,           0.44]);
fprintf('Share of L of top 0.05        = %9.2f  %9.2f \n',     [fRtopL5,           0.64]);
fprintf('Share of L of top 0.1         = %9.2f  %9.2f \n',     [fRtopL10,          0.73]);
fprintf('Share of L of top 0.2         = %9.2f  %9.2f \n',     [fRtopL20,          0.82]);
fprintf('\n');
fprintf('Share of K of top 0.01        = %9.2f  %9.2f \n',     [fRtopK1,           0.62]);
fprintf('Share of K of top 0.05        = %9.2f  %9.2f \n',     [fRtopK5,           0.78]);
fprintf('Share of K of top 0.1         = %9.2f  %9.2f \n',     [fRtopK10,          0.85]);
fprintf('Share of K of top 0.2         = %9.2f  %9.2f \n',     [fRtopK20,          0.91]);

fprintf('\n');
fprintf('\n');
fprintf('Share of Y of top 0.01 young  = %9.2f  %9.2f \n',       [fRtopYy1,            0.34]);
fprintf('Share of Y of top 0.05 young  = %9.2f  %9.2f \n',       [fRtopYy5,            0.53]);
fprintf('Share of Y of top 0.1  young  = %9.2f  %9.2f \n',       [fRtopYy10,           0.64]);
fprintf('Share of Y of top 0.2  young  = %9.2f  %9.2f \n',       [fRtopYy20,           0.75]);
fprintf('\n');

fprintf('Share of Y of top 0.01 middle = %9.2f  %9.2f \n',       [fRtopYm1,            0.42]);
fprintf('Share of Y of top 0.05 middle = %9.2f  %9.2f \n',       [fRtopYm5,            0.63]);
fprintf('Share of Y of top 0.1  middle = %9.2f  %9.2f \n',       [fRtopYm10,           0.73]);
fprintf('Share of Y of top 0.2  middle = %9.2f  %9.2f \n',       [fRtopYm20,           0.83]);
fprintf('\n');

fprintf('Share of Y of top 0.01  old   = %9.2f  %9.2f \n',       [fRtopYo1,            0.53]);
fprintf('Share of Y of top 0.05  old   = %9.2f  %9.2f \n',       [fRtopYo5,            0.77]);
fprintf('Share of Y of top 0.1   old   = %9.2f  %9.2f \n',       [fRtopYo10,           0.86]);
fprintf('Share of Y of top 0.2   old   = %9.2f  %9.2f \n',       [fRtopYo20,           0.92]);

fprintf('\n');
fprintf('\n');

fprintf('Share of L of top 0.01 young  = %9.2f  %9.2f \n',       [fRtopLy1,            0.23]);
fprintf('Share of L of top 0.05 young  = %9.2f  %9.2f \n',       [fRtopLy5,            0.41]);
fprintf('Share of L of top 0.1  young  = %9.2f  %9.2f \n',       [fRtopLy10,           0.53]);
fprintf('Share of L of top 0.2  young  = %9.2f  %9.2f \n',       [fRtopLy20,           0.67]);
fprintf('\n');

fprintf('Share of L of top 0.01 middle = %9.2f  %9.2f \n',       [fRtopLm1,            0.31]);
fprintf('Share of L of top 0.05 middle = %9.2f  %9.2f \n',       [fRtopLm5,            0.52]);
fprintf('Share of L of top 0.1  middle = %9.2f  %9.2f \n',       [fRtopLm10,           0.63]);
fprintf('Share of L of top 0.2  middle = %9.2f  %9.2f \n',       [fRtopLm20,           0.76]);
fprintf('\n');

fprintf('Share of L of top 0.01  old   = %9.2f  %9.2f \n',       [fRtopLo1,            0.45]);
fprintf('Share of L of top 0.05  old   = %9.2f  %9.2f \n',       [fRtopLo5,            0.70]);
fprintf('Share of L of top 0.1   old   = %9.2f  %9.2f \n',       [fRtopLo10,           0.80]);
fprintf('Share of L of top 0.2   old   = %9.2f  %9.2f \n',       [fRtopLo20,           0.89]);
fprintf('\n');
fprintf('\n');

fprintf('Share of K of top 0.01 young  = %9.2f  %9.2f \n',       [fRtopKy1,            0.50]);
fprintf('Share of K of top 0.05 young  = %9.2f  %9.2f \n',       [fRtopKy5,            0.66]);
fprintf('Share of K of top 0.1  young  = %9.2f  %9.2f \n',       [fRtopKy10,           0.75]);
fprintf('Share of K of top 0.2  young  = %9.2f  %9.2f \n',       [fRtopKy20,           0.84]);
fprintf('\n');
fprintf('Share of K of top 0.01 middle = %9.2f  %9.2f \n',       [fRtopKm1,            0.56]);
fprintf('Share of K of top 0.05 middle = %9.2f  %9.2f \n',       [fRtopKm5,            0.72]);
fprintf('Share of K of top 0.1  middle = %9.2f  %9.2f \n',       [fRtopKm10,           0.80]);
fprintf('Share of K of top 0.2  middle = %9.2f  %9.2f \n',       [fRtopKm20,           0.87]);
fprintf('\n');
fprintf('Share of K of top 0.01  old   = %9.2f  %9.2f \n',       [fRtopKo1,            0.59]);
fprintf('Share of K of top 0.05  old   = %9.2f  %9.2f \n',       [fRtopKo5,            0.80]);
fprintf('Share of K of top 0.1   old   = %9.2f  %9.2f \n',       [fRtopKo10,           0.88]);
fprintf('Share of K of top 0.2   old   = %9.2f  %9.2f \n',       [fRtopKo20,           0.93]);





