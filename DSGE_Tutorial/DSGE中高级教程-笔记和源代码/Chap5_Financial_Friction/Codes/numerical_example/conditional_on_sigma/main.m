clear
close all
dbstop if error

R=1.011470654984364;
Rk=1.018806396097570;
sp=Rk/R;%risk spread
mu=0.2148945111281970;
% Fomegabar=0.0055885692972290; %bankruptancy rate are set. 这个其实是用不到的
sigma =  0.259192677224171*3;

%when given sigma, mu and spread, find an omega that solves the efficiency
%condition, i.e., the foc of the entrepreneurial utility maximization
%problem.

[omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_sigma(sigma,mu,sp);

leverage=1/(1-Gam_muG*sp);%bank zero profit conditon
spread=sp*omega*leverage/(leverage-1);%interest rate spread
spreadAR=400*(spread-1); %annulized rate
util=(1-Gamma)*sp*leverage; %utility

%% graph the parameters of the equilibrium contract for various settings of
%the parameters....
%1. mu  2.risk spread  3. sigma
lev=[];
Pe=[]; % risk spread
%now, compute the equilibrium starting from just the primitive parameters
mmu=.02:.005:.33;
for ii = 1:length(mmu)
    xmu=mmu(ii);
    [omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_sigma(sigma,xmu,sp);
    lev(ii)=1/(1-Gam_muG*sp); %bank zero profit condition
    Pe(ii)=sp*omega*lev(ii)/(lev(ii)-1);   %interest rate spread  
    bank(ii)=F; %bankruptancy rate
end
ia=3;
ib=3;
subplot(ia,ib,1)
plot(mmu,lev)
ylabel('leverage')
xlabel('\mu')
axis tight

subplot(ia,ib,2)
plot(mmu,400*(Pe-1))
ylabel('400*(Z/R-1)')
xlabel('\mu')
axis tight

subplot(ia,ib,3)
plot(mmu,100*bank)
ylabel('bankruptcy rate (%)')
xlabel('\mu')
axis tight

bank=[];
lev=[];
Pe=[];
ssp1=1.001:.0025:1.03;
for ii = 1:length(ssp1)
    xsp1=ssp1(ii);
    [omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_sigma(sigma,mu,xsp1);
    lev(ii)=1/(1-Gam_muG*xsp1);
    Pe(ii)=sp*omega*lev(ii)/(lev(ii)-1);
    bank(ii)=F;
end
subplot(ia,ib,4)
plot(400*(ssp1-1),lev)
ylabel('leverage')
xlabel('400*(Rk/R-1)')
axis tight

subplot(ia,ib,5)
plot(400*(ssp1-1),400*(Pe-1))
ylabel('400*(Z/R-1)')
xlabel('400*(Rk/R-1)')
axis tight

subplot(ia,ib,6)
plot(400*(ssp1-1),100*bank)
ylabel('bankruptcy rate (%)')
xlabel('400*(Rk/R-1)')
axis tight

bank=[];
lev=[];
Pe=[];
ssigma=0.03:.01:1;
for ii = 1:length(ssigma)
    xsigma=ssigma(ii);
    [omega,G,F,Gamma,Gam_muG,Fprime,er] = get_omega_cond_sigma(xsigma,mu,sp);
    lev(ii)=1/(1-Gam_muG*sp);
    Pe(ii)=sp*omega*lev(ii)/(lev(ii)-1);
    bank(ii)=F;
end
subplot(ia,ib,7)
plot(ssigma,lev)
ylabel('leverage')
xlabel('\sigma')
axis tight

subplot(ia,ib,8)
plot(ssigma,400*(Pe-1))
ylabel('400*(Z/R-1)')
xlabel('\sigma')
axis tight

subplot(ia,ib,9)
plot(ssigma,100*bank)
ylabel('bankruptcy rate (%)')
xlabel('\sigma')
axis tight

text='Leverage and interest rate spread, for alternative parameter settings';
sgtitle(text);