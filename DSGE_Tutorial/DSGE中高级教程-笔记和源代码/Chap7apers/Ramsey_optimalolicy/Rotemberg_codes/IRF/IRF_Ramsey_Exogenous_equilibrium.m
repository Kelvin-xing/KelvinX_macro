%this file will plot the RBC model's IRF from exogenous monetary policy
%rule and Ramsey equilibrium.

clear all;
close all;
clc;

%we iterate on parameter nu
epsilon = 5;
nu_arr = [0, 1/(epsilon-1)];
%run the dynare mod file and store the results;
for ii=1:length(nu_arr)
    nu = nu_arr(ii);
    save parameterfile_irf nu;    
    dynare exogenous_monetary_rule noclearall;
    str=['save exogenous_mon_rule',int2str(ii)];
    eval(str);
  
    dynare Rotemberg_ramsey_policy noclearall;
    str=['save ramsey_policy',int2str(ii)];
    eval(str);
end
%% plot the IRF together
close all;
clear all;

epsilon = 5;
nu_arr = [0, 1/(epsilon-1)];

irfn=40;
T=1:1:irfn;
netinflation=zeros(length(nu_arr),irfn);
hours = zeros(length(nu_arr),irfn);
net_nom_interest =zeros(length(nu_arr),irfn);
consumption = zeros(length(nu_arr),irfn);
technology = zeros(length(nu_arr),irfn);
real_rate =zeros(length(nu_arr),irfn);


for ii=1:length(nu_arr)
    str=['load exogenous_mon_rule',int2str(ii)];
    eval(str);
    netinflation(1,:)=400*(pie_eps_A + pie_SS-1);
    hours(1,:) = 100*h_eps_A/h_SS;
    net_nom_interest(1,:) = 400*(R_eps_A +R_SS);
    consumption(1,:) = 100*C_eps_A/C_SS;
    technology(1,:) = A_eps_A;
    real_rate(1,:) = 400*(R_eps_A +R_SS- pie_eps_A-pie_SS+1);
    
    str=['load ramsey_policy',int2str(ii)];
    eval(str);
    netinflation(2,:)=400*(pie_eps_A + pie_SS-1);
    hours(2,:) = 100*h_eps_A/h_SS;
    net_nom_interest(2,:) = 400*(R_eps_A +R_SS);
    consumption(2,:) = 100*C_eps_A/C_SS;
    technology(2,:) = A_eps_A;
    real_rate(2,:) = 400*(R_eps_A+R_SS- pie_eps_A -pie_SS +1);
    
    figure;
    subplot(321);
    plot(T,netinflation(1,:),'k*-',T,netinflation(2,:));
    legend('exog. mon. rule','Ramsey')
    ylabel('percent')
    title('net inflation (APR)')
    axis tight
    
    subplot(322);
    plot(T,hours(1,:),'k*-',T,hours(2,:));
    ylabel('percent of steady state')
    title('hours worked')
    axis tight
    
    subplot(323);
    plot(T,net_nom_interest(1,:),'k*-',T,net_nom_interest(2,:));
    ylabel('percent')
    title('net nominal interest rate (APR)')
    axis tight
    
    subplot(324);
    plot(T,consumption(1,:),'k*-',T,consumption(2,:));
    ylabel('percent of steady state')
    title('consumption')
    axis tight
    
    subplot(325);
    plot(T,technology(1,:),'k*-',T,technology(2,:));
    ylabel('unit of standard deviation')
    title('technology')
    axis tight
    
    subplot(326);
    plot(T,real_rate(1,:),'k*-',T,real_rate(2,:));
    ylabel('percent')
    title('deviation of real rate from s.s.(APR)')
    axis tight    
end