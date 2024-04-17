%price dispersion and trend inflation rate
clear;
pis = 0.98:0.001:1.02;
phi =.75;
epsilon =10;
for ii=1:length(pis)
    pisharps(ii) = ((pis(ii)^(1-epsilon) - phi)/(1-phi))^(1/(1-epsilon));
    vps(ii) = (1-phi)*(pis(ii)/pisharps(ii))^epsilon/(1- pis(ii)^epsilon*phi);
end

plot(pis,vps,'LineWidth',2)
axis tight
ylabel('v^{p}_{t}')
xlabel('\pi')
title('Steady state price dispersion and trend inflation')