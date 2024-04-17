%% spurious steady state, 2015-10-27@BJ
%Written by Xiangyang Li.
%2nd approximation under gamma=2
ks = 3.8774;
kt = -0.2: 0.01:3.5; %deviation from steady states;
yt = zeros(length(kt),1)';
xt = zeros(length(kt),1)';
yt(1) = kt(1);
gk = 0.98;
ga = 0.063;
gkk = 0.014;
gaa= 0.067;
gss = 0.000024;
gka = -0.035;
Vara = 0.0320^2; %var(a) = Veps/(1-rho^2); rho = 0.95; Veps =0.01^2;
at = [0,2*sqrt(Vara),4*sqrt(Vara)];
for jj=1:length(at)
    for ii=2:length(kt)
        yt(ii)= gk*kt(ii-1)+ga*at(jj)+...
        0.5*(gkk*kt(ii-1)^2+gaa*at(jj)^2+gss)+gka*kt(ii-1)*at(jj);  
        xt(ii)=yt(ii) - kt(ii-1);
    end
    switch jj
        case 1
            plot(kt(2:end-1),xt(3:end));
        case 2
            plot(kt(2:end-1),xt(3:end),'r');
        case 3
            plot(kt(2:end-1),xt(3:end),'g');
    end 
    hold on;
end
plot(kt(2:end-1),zeros(length(kt(2:end-1))),'cyan');
hold off;
axis tight;