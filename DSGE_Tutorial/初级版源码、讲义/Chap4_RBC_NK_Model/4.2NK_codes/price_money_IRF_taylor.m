%% Sticky price Taylor rule for technology shock
% back out the nominal Price level and nominal money
% this scripts should be run immediately after running mod file
T = 40;
lnP_ea(1) = pi_ea(1);
M_ea(1) = m_ea(1)+lnP_ea(1);
for ii=2:T
    lnP_ea(ii) = lnP_ea(ii-1) + pi_ea(ii); %in log, intial price level =1
    M_ea(ii) = m_ea(ii)+lnP_ea(ii); %in log
end
figure;
subplot(2,2,1);
plot(1:1:T,lnP_ea);
title('Price level');
subplot(2,2,2);
plot(1:1:T,M_ea);
title('Nominal Money');
subplot(2,2,3);
plot(1:1:T,outputgap_ea);
title('Output Gap');
subplot(2,2,4);
plot(1:1:T,a_ea);
title('technology');
%% sticky price, taylor rule, for monetary policy shock
T = 40;
lnP_ei(1) = pi_ei(1);
M_ei(1) = m_ei(1)+lnP_ei(1);
for ii=2:T
    lnP_ei(ii) = lnP_ei(ii-1) + pi_ei(ii); %in log, intial price level =1
    M_ei(ii) = m_ei(ii)+lnP_ei(ii); %in log
end
figure;
subplot(2,2,1);
plot(1:1:T,lnP_ei);
title('Price level ');
subplot(2,2,2);
plot(1:1:T,M_ei);
title('Nominal Money');
subplot(2,2,3);
plot(1:1:T,yf_ei);
axis([0 40 -0.1 0.1])
title('Flexible output');