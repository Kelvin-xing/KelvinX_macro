Tinit = 50;       % periods to initialize distribution of assets
T     = 6;        % periods to simulate allocations
Nf    = 50000;    % # producers

randn('state', 100);
rand('state',  100);

eesim    = zeros(Nf,T);       % transitory component
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

% run a bit without saving data to get initial distribution

for t = 1 : Tinit+T

    if t==1
         
        e      = esim(:,t);
        
        state  = [smin(1)*ones(Nf,1), e]; 
        a      = state(:,1);
        x      = funeval(cx,fspace,state);

    else
       
        e      = esim(:,t);                 

        state  = [x, e];
        
        a      = state(:,1);
        x      = max(funeval(cx, fspace, state), smin(1));
            
          if t > Tinit
               
              asim(:,t-Tinit)    = a;
              eesim(:,t-Tinit)   = e;
              xsim(:,t-Tinit)    = x;    
              
          end
    end
end

a    = asim';
x    = xsim';
e    = eesim';
rr   = reshape(funeval(cr,fspace,[a(:), e(:)]), size(a,1), size(a,2));
e    = egrid(e);


L   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(rr+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e);
K   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(rr+delta).^((alpha*eta-1)/(1-eta)).*exp(e);
Y   = exp(e).^(1-eta).*(L.^alpha.*K.^(1-alpha)).^eta;
D   = Y - W.*L - (r+delta).*K;  

cons = D + (1+r)*a - x;                    % consumption, (recall x is total savings)

% Frictionless allocations

Lf   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*(r+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e);
Kf   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*(r+delta).^((alpha*eta-1)/(1-eta)).*exp(e);
Yf   = exp(e).^(1-eta).*(Lf.^alpha.*Kf.^(1-alpha)).^eta;

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

DtoYm = sum(Debt(:))./sum(Y(:));

fprintf('aggregate Debt to GDP         = %9.2f  %9.2f \n',       [DtoYm,       DtoY]); 
%fprintf('aggregate posit. Debt to GDP  = %9.2f  %9.2f \n',       [sum(Debt(Debt>0))./sum(Y(:)),       DtoY]); 

%fprintf('\n');
%fprintf('Variance exog component       = %9.2f \n',       varz);

