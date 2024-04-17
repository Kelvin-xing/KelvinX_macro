%2015-10-14@Beijing
%written by Xiangyang  Li
%% extended example for  Prunning -1 
clear;
rho =0.8;
alpha =.5;
sigma = 0.1;
N=500; %number of simulation
z=-0.6:0.01:1;
eps = zeros(length(z),1);
for ii=1:(length(z))
    y(ii)= rho*z(ii)+alpha*z(ii)^2+eps(ii);
    x(ii)=rho*z(ii)+alpha*z(ii)^2 + 2*sigma;
    w(ii)=rho*z(ii)+alpha*z(ii)^2 - 2*sigma;
end
plot(z,y,'g');
hold on;
plot(z,z,'cyan','LineWidth',2);
hold on;
plot(z,x,':r');
hold on;
plot(z,w,'b.')

%% Pruning - for small shock
eps = sigma/10*randn(N,1); 
y(1)=0;
x(1)=0;
for jj=2:N
    y(jj)= rho*y(jj-1)+alpha*y(jj-1)^2+eps(jj);
    x(jj)= rho*x(jj-1)+eps(jj);
end
% the two line alomst the same for small variance but differ tremendously
% when variance is large.
plot(1:N,y,'-r.',1:N,x)

%the difference btw y and x;
dist=y-x;
%% Pruning - for large shock
eps = sigma*randn(N,1); 
y(1)=0;
x(1)=0;
for jj=2:N
    y(jj)= rho*y(jj-1)+alpha*y(jj-1)^2+eps(jj);
    x(jj)= rho*x(jj-1)+eps(jj);
end
% the two line alomst the same for small variance but differ tremendously
% when variance is large.
plot(1:N,y,'-r.',1:N,x)
ylim([-0.5 2])
%% spurious steady state 
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
        yt(ii)= gk*kt(ii-1)+ga*at(jj)+0.5*(gkk*kt(ii-1)^2+gaa*at(jj)^2+gss)+gka*kt(ii-1)*at(jj);  
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
%% the root of spurious steady states
clear;
clc;
gk = 0.98;
ga = 0.063;
gkk = 0.014;
gaa= 0.067;
gss = 0.000024;
gka = -0.035;
Vara = 0.0320^2;
sigma =1 ;
kstar = 3.87;
at = [0,2*sqrt(Vara),4*sqrt(Vara)];
for jj=1:length(at)
    ca = 0.5*gkk;
    cb = gk+gka*at(jj) -1;
    cc = ga*at(jj)+0.5*(gaa*at(jj)^2+gss);
    g1 = (-cb + sqrt(cb^2 - 4*ca*cc))/2/ca
    g2 = (-cb - sqrt(cb^2 - 4*ca*cc))/2/ca
    K1 = exp(kstar +g1)
    K2 = exp(kstar +g2)
end
%% price distortion under sticky price
%which are asymptotical to unity.
clear;
theta = 0.75;
eps = 10;
p1=zeros(30,1);
%you can try different initial values;
p1(1) = 0.8;
p2=zeros(30,1);
p2(1) = 0.8;
for ii=2:length(p1)
    p2(ii)= (1-theta+theta*p2(ii-1)^(1-eps))^(1/(1-eps));
    
    %the same as above, but converge slower;
    p1(ii)= (1-theta+theta*p1(ii-1)^(eps-1))^(1/(eps-1));
end
%converge to unity no matter p(1) is less or larger than unity;
plot(1:30,p1,'-r*',1:30,p2,':b.');







