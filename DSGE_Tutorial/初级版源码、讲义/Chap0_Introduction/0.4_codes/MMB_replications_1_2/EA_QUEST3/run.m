% Run and plot IRF of replication file 

% EA_QUEST3 model

clear all;
clc;

% Adjust path to folder where replication file is stored
cd([cd '/EA_QUEST3_rep']);

% Run replication dynare file
dynare EA_Quest3_rep


infirf = E_PHIC_E_EPS_M;
intirf = E_INOM_E_EPS_M;  
yirf =  E_LYGAP_E_EPS_M; 
infairf = E_PHI_E_EPS_M;


% Go back to original path
cd('..');


% Plot replicated IRF
t = 1:1:length(infirf);
zeroline = ones(length(t),1)*0;

subplot(2,2,1);
plot(t,yirf,'LineWidth',2);
axis([0 42 -3*10^(-3) 1*10^(-3)]);
xlabel('quarters','FontSize',8);
title('Output gap: log(YGAP_{t})','FontSize',10);

subplot(2,2,2);
plot(t,infairf,'LineWidth',2);
axis([0 42 -8*10^(-4) 0]);
xlabel('quarters','FontSize',8);
title('Annual inflation: \pi_{t}','FontSize',10);

subplot(2,2,3);
plot(t,intirf,'LineWidth',2);
axis([0 42 -5*10^(-4) 15*10^(-4)]);
xlabel('quarters','FontSize',8);
title('Interest Rate: inom_{t}','FontSize',10);

subplot(2,2,4);
plot(t,infirf,'LineWidth',2);
axis([0 42 -1*10^(-3) 0]);
xlabel('quarters','FontSize',8);
%legend1 = legend('replication','original');
%set(legend1,'Position',[0.6435 0.2317 0.2089 0.1]);
title('Quarterly inflation: \pi_{t}^{c}','FontSize',10);




