Ys=[Ys1 Ys2 Ys3];
Y=Ys;

Ls=[Ls1 Ls2 Ls3];
L=Ls;

Ks=[Ks1 Ks2 Ks3];
K=Ks;

es=[es1 es2 es3];
e=es;

rs=[rs1 rs2 rs3];
rr=rs;


% wedgel = (rr+delta).^(-(1-alpha)*eta/(1-eta)).*mean(exp(e(:)))./mean((rr(:)+delta).^(-(1-alpha)*eta/(1-eta)).*exp(e(:)));
% wedgek = (rr+delta).^((alpha*eta-1)/(1-eta)).*mean(exp(e(:)))./mean((rr(:)+delta).^((alpha*eta-1)/(1-eta)).*exp(e(:)));
% 
% TFPratio = mean(wedgel(:).^(alpha*eta).*wedgek(:).^((1-alpha)*eta).*exp(e(:)))./mean(exp(e(:)));
% 
% % autocorrelation Y for firms in modern sector
% 
Yc   = vec(Y(2:end,:));
Yp   = vec(Y(1:end-1,:));
% Lc   = vec(L(2:end,:));
% Lp   = vec(L(1:end-1,:));
% Kc   = vec(K(2:end,:));
% Kp   = vec(K(1:end-1,:));
% 
% 
acorY = corr([log(Yc),log(Yp)]);
% acorL = corr([log(Lc),log(Lp)]);
% acorK = corr([log(Kc),log(Kp)]);
% 
Yc = vec(Y(3:end,:));
Yp = vec(Y(1:end-2,:));
% Lc = vec(L(3:end,:));
% Lp = vec(L(1:end-2,:));
% Kc = vec(K(3:end,:));
% Kp = vec(K(1:end-2,:));
% 
%acorY2 = corr([log(Yc), log(Yp)]);
% acorL2 = corr([log(Lc), log(Lp)]);
% acorK2 = corr([log(Kc), log(Kp)]);
% 
Yc = vec(Y(4:end,:));
Yp = vec(Y(1:end-3,:));
% Lc = vec(L(4:end,:));
% Lp = vec(L(1:end-3,:));
% Kc = vec(K(4:end,:));
% Kp = vec(K(1:end-3,:));
% 
% 
acorY3 = corr([log(Yc), log(Yp)]);
% acorL3 = corr([log(Lc), log(Lp)]);
% acorK3 = corr([log(Kc), log(Kp)]);
% 
% Yc = vec(Y(5:end,:));
% Yp = vec(Y(1:end-4,:));
% Lc = vec(L(5:end,:));
% Lp = vec(L(1:end-4,:));
% Kc = vec(K(5:end,:));
% Kp = vec(K(1:end-4,:));
% 
% acorY4 = corr([log(Yc), log(Yp)]);
% acorL4 = corr([log(Lc), log(Lp)]);
% acorK4 = corr([log(Kc), log(Kp)]);
% 
Yc = vec(Y(6:end,:));
Yp = vec(Y(1:end-5,:));
% Lc = vec(L(6:end,:));
% Lp = vec(L(1:end-5,:));
% Kc = vec(K(6:end,:));
% Kp = vec(K(1:end-5,:));
% 
acorY5 = corr([log(Yc), log(Yp)]);
% acorL5 = corr([log(Lc), log(Lp)]);
% acorK5 = corr([log(Kc), log(Kp)]);
% 
% % growth rates:
% 
Yc = vec(Y(2:end,:));
Yp = vec(Y(1:end-1,:));
% Lc = vec(L(2:end,:));
% Lp = vec(L(1:end-1,:));
% Kc = vec(K(2:end,:));
% Kp = vec(K(1:end-1,:));
% 
dY  = vec(log(Yc./Yp)); 
% dL  = vec(log(Lc./Lp));
% dK  = vec(log(Kc./Kp));
% 
% 
% 
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
stdY=std(dY(:));
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
fprintf('s.d.  dlog(Y)                 = %9.2f  %9.2f \n',      [std(dY(:)),           sdYdata]);
fprintf('s.d.   log(Y)                 = %9.2f  %9.2f \n',       [sY,                   sYdata]);
fprintf('\n');
fprintf('cor y y1                      = %9.2f  %9.2f \n',       [ac1,                  ac1Ydata]);
fprintf('cor y y3                      = %9.2f  %9.2f \n',       [ac3,                  ac3Ydata]);
fprintf('cor y y5                      = %9.2f  %9.2f \n',       [ac5,                  ac5Ydata]);

fprintf('\n');
% 
% DtoYm = sum(Debt(:))./sum(Y(:));
% 
% fprintf('aggregate Debt to GDP         = %9.2f  %9.2f \n',       [DtoYm,       DtoY]); 
%fprintf('aggregate posit. Debt to GDP  = %9.2f  %9.2f \n',       [sum(Debt(Debt>0))./sum(Y(:)),       DtoY]); 

%fprintf('\n');
%fprintf('Variance exog component       = %9.2f \n',       varz);


alpha = 1/3*(alpha1 + alpha2 + alpha3);

TFPratio = mean(exp(e(:)).*(Y(:)./K(:)).^(-(1-alpha)*eta/(1-eta))).^(1-alpha*eta)/...
    mean(exp(e(:)).*(Y(:)./K(:)).^((alpha*eta-1)/(1-eta))).^((1-alpha)*eta)/mean(exp(e(:))).^(1-eta);


fprintf('Naive Misallocation Loss Modern     = %9.3f \n',          -log(TFPratio)*100);

