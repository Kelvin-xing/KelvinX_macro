Tinit = 50;       % periods to initialize distribution of assets
T     = 6;        % periods to simulate allocations
Nf    = 50000;    % # producers

randn('state', 100);
rand('state',  100);

eesim    = zeros(Nf,T);       % transitory component
asim     = zeros(Nf,T);       % beginning of period assets
xsim     = zeros(Nf,T);        
rrsim    = zeros(Nf,T);
E1sim    = zeros(Nf,T);

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

% run a bit without saving data to get initial distribution

for t = 1 : Tinit+T

    if t==1
         
        e      = esim(:,t);
        
        state  = [smin(1)*ones(Nf,1), e, zeros(Nf,1)]; 
        a      = state(:,1);
        x      = funeval(cx,fspace,state);

    else
       
        e      = esim(:,t);                 

        state  = [x, e, esim(:,t-1)];
        
        a      = state(:,1);
        x      = max(funeval(cx, fspace, state), smin(1));
            
          if t > Tinit
               
              asim(:,t-Tinit)    = a;
              eesim(:,t-Tinit)   = e;
              xsim(:,t-Tinit)    = x;
              
              Phi = funbas(fspace,state); 
              
              rrsim(:,t-Tinit)   = Phi*cr;
              E1sim(:,t-Tinit)    = Phi*ce;
              
              
          end
    end
end

a    = asim';
x    = xsim';
e    = eesim';
rr   = rrsim';
E1   = E1sim';

e    = egrid(e);

K  = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*E1.^((1-alpha*eta)/(1-eta));
Kf = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*E1.^((1-alpha*eta)/(1-eta));

L  = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e).^((1-eta)./(1-alpha*eta)).*E1.^((1-alpha)*eta/(1-eta));
Lf = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e).^((1-eta)./(1-alpha*eta)).*E1.^((1-alpha)*eta/(1-eta));

Y  = exp(e).^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
Yf = exp(e).^(1-eta).*(Lf.^alpha.*Kf.^(1-alpha)).^eta;

D = Y - W.*L - (r+delta).*K; 

cons = D + (1+r)*a - x;                    % consumption, (recall x is total savings)

A = a;

Debt   = K - A;

% three ways to compute TFP losses in Modern Sector (TFP calculation only
% up to a scalar of normalization, may differ from number in ergodic.m)

TFP      = mean(Y(:))./(mean(L(:)).^alpha*mean(K(:)).^(1-alpha))^eta;
TFPbest  = mean(Yf(:))./(mean(Lf(:)).^alpha*mean(Kf(:)).^(1-alpha))^eta;
TFP2best = mean(exp(e(:))).^(1-eta);

wedgel = (rr+delta).^(-(1-alpha)*eta/(1-eta)).*mean(exp(e(:)))./mean((rr(:)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e(:)));
wedgek = (rr+delta).^((alpha*eta-1)/(1-eta)).*mean(exp(e(:)))./mean((rr(:)+delta).^((alpha*eta-1)/(1-eta)).*exp(e(:)));

TFPratio = mean(wedgel(:).^(alpha*eta).*wedgek(:).^((1-alpha)*eta).*exp(e(:)))./mean(exp(e(:)));

% autocorrelation Y for firms in modern sector

Yc   = vec(Y(2:end,:));
Yp   = vec(Y(1:end-1,:));
Lc   = vec(L(2:end,:));
Lp   = vec(L(1:end-1,:));
Kc   = vec(K(2:end,:));
Kp   = vec(K(1:end-1,:));


acorY = corr([log(Yc),log(Yp)]);
acorL = corr([log(Lc),log(Lp)]);
acorK = corr([log(Kc),log(Kp)]);

Yc = vec(Y(3:end,:));
Yp = vec(Y(1:end-2,:));
Lc = vec(L(3:end,:));
Lp = vec(L(1:end-2,:));
Kc = vec(K(3:end,:));
Kp = vec(K(1:end-2,:));

acorY2 = corr([log(Yc), log(Yp)]);
acorL2 = corr([log(Lc), log(Lp)]);
acorK2 = corr([log(Kc), log(Kp)]);

Yc = vec(Y(4:end,:));
Yp = vec(Y(1:end-3,:));
Lc = vec(L(4:end,:));
Lp = vec(L(1:end-3,:));
Kc = vec(K(4:end,:));
Kp = vec(K(1:end-3,:));


acorY3 = corr([log(Yc), log(Yp)]);
acorL3 = corr([log(Lc), log(Lp)]);
acorK3 = corr([log(Kc), log(Kp)]);

Yc = vec(Y(5:end,:));
Yp = vec(Y(1:end-4,:));
Lc = vec(L(5:end,:));
Lp = vec(L(1:end-4,:));
Kc = vec(K(5:end,:));
Kp = vec(K(1:end-4,:));

acorY4 = corr([log(Yc), log(Yp)]);
acorL4 = corr([log(Lc), log(Lp)]);
acorK4 = corr([log(Kc), log(Kp)]);

Yc = vec(Y(6:end,:));
Yp = vec(Y(1:end-5,:));
Lc = vec(L(6:end,:));
Lp = vec(L(1:end-5,:));
Kc = vec(K(6:end,:));
Kp = vec(K(1:end-5,:));

acorY5 = corr([log(Yc), log(Yp)]);
acorL5 = corr([log(Lc), log(Lp)]);
acorK5 = corr([log(Kc), log(Kp)]);

% growth rates:

Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));
Lc = vec(L(2:end,:));
Lp = vec(L(1:end-1,:));
Kc = vec(K(2:end,:));
Kp = vec(K(1:end-1,:));

dY  = vec(log(Yc./Yp)); 
dL  = vec(log(Lc./Lp));
dK  = vec(log(Kc./Kp));



% Pin down variance of z needed to match the size distributions

if findvarz

varo  = var(e(:));
varz  = (0.0: .01: 2)';
varY  = var(log(Y(:)));
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
    
varo  = var(e(:));
varY  = var(log(Y(:)));
sdY   = sqrt(varY + varz);
ac1   = (acorY.*varY  + varz)./(varY + varz);
ac3   = (acorY3.*varY + varz)./(varY + varz);
ac5   = (acorY5.*varY + varz)./(varY + varz);
    
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
fprintf('aggregate Debt to GDP         = %9.2f  %9.2f \n',       [sum(Debt(:))./sum(Y(:)),       DtoY]); 
fprintf('aggregate posit. Debt to GDP  = %9.2f  %9.2f \n',       [sum(Debt(Debt>0))./sum(Y(:)),       DtoY]); 

fprintf('\n');
fprintf('Variance exog component       = %9.2f \n',       varz);



fprintf('\n');
fprintf('\n');
fprintf('\n');

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(K(:)),   mean(A(:)) + Aw]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [mean(L(:)),   Lbar]);
fprintf('\n');

fprintf('\n');
disp('Model Implications')
fprintf('\n');

Da = mean(K(:)) - mean(A(:)) - Aw;       % Aggregate debt position

fprintf('Output                        = %9.3f \n',          mean(Y(:))/mean(L(:)));
fprintf('Consumption                   = %9.3f \n',          (mean(Y(:)) - delta*mean(K(:)) - r*Da)/mean(L(:)));
fprintf('Investment                    = %9.3f \n',          (delta*mean(K(:)))/mean(L(:)));
fprintf('TFP Modern                    = %9.3f \n',          sum(Y(:))./(sum(L(:)).^alpha*sum(K(:)).^(1-alpha))^eta*numel(Y)^(eta-1)); % express in per capita units, comparable to ergodic
fprintf('Misallocation Loss Modern     = %9.3f \n',          log(TFPbest/TFP)*100);


fprintf('\n');


ract = rr(:);

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


varL  = var(log(L(:)));
sdL   = sqrt(varL + varz);
acL1  = (acorL.*varL  + varz)./(varL + varz);
acL3  = (acorL3.*varL + varz)./(varL + varz);
acL5  = (acorL5.*varL + varz)./(varL + varz);

varK  = var(log(K(:)));
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


% Elasticity of y/k to y

y = log(Y); 
k = log(K);

yk = y - k;

dy  = vec(y(2:end,:) - y(1:end-1,:));
cyk = vec(yk(2:end,:));
lyk = vec(yk(1:end-1,:));


xx = [ones(size(cyk)), dy, lyk];
yy = cyk;

beta = inv(xx'*xx)*xx'*yy;
fprintf('Elasticity of yk to dy given yk(-1) = %9.2f  %9.2f \n',     [beta(2),         0.72]);



de = vec(e(2:end,:)-e(1:end-1,:));
xx = [ones(size(cyk)), de, lyk];
yy = cyk;

beta = inv(xx'*xx)*xx'*yy;
fprintf('Elasticity of yk to de given yk(-1) = %9.2f  %9.2f \n',     [beta(2),         0.72]);
fprintf('\n');
fprintf('S.d. ln(Y/K)  = %9.2f  \n',     [std(yk(:))]);
fprintf('\n');

fprintf('\n');
TFPratio = mean(exp(e(:)).*(Y(:)./K(:)).^(-(1-alpha)*eta/(1-eta))).^(1-alpha*eta)/...
    mean(exp(e(:)).*(Y(:)./K(:)).^((alpha*eta-1)/(1-eta))).^((1-alpha)*eta)/mean(exp(e(:))).^(1-eta);


fprintf('Naive Misallocation Loss Modern     = %9.3f \n',          -log(TFPratio)*100);

% Some implications that can test with data

fprintf('\n');
fprintf('\n');

fprintf('s.d.  dlog(Y)                 = %9.2f \n',     std(dY(:)));
fprintf('s.d.  dlog(K)                 = %9.2f \n',     std(dK(:)));
fprintf('\n');

Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));

Kc = vec(K(2:end,:));
Kp = vec(K(1:end-1,:));

ec = vec(e(2:end,:));
ep = vec(e(1:end-1,:));

dY  = log(Yc(:)./Yp(:)); 
dK  = log(Kc(:)./Kp(:));
de  = ec(:) - ep(:);


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


YK = log(Y(:)./K(:)); 
e  = e(:); 

fprintf('var  YK                       = %9.3f \n',    var(YK));


TFPloss = mean(exp(e - (1-alpha)*eta/(1-eta)*YK))^(1-alpha*eta)/...
          mean(exp(e + (alpha*eta-1)/(1-eta)*YK))^((1-alpha)*eta)/mean(exp(e))^(1-eta);

      
fprintf('\n');

fprintf('TFP loss                      = %9.3f \n',    -log(TFPloss)*100);

