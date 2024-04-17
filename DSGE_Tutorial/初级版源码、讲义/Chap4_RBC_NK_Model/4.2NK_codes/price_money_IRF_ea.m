%% back out the nominal Price level and nominal money to technology shock
% this scripts should be run immediately after running mod file
T = 40;
lnP_em(1) = pi_ea(1);
M_em(1) = m_ea(1)+lnP_em(1);
for ii=2:T
    lnP_em(ii) = lnP_em(ii-1) + pi_ea(ii); %in log, intial price level =1
    M_em(ii) = m_ea(ii)+lnP_em(ii); %in log
end
figure;
subplot(1,3,1);
plot(1:1:T,lnP_em);
title('Price level ');
subplot(1,3,2);
plot(1:1:T,M_em);
axis([0 40 -0.1 0.1])
title('Nominal Money');
subplot(1,3,3);
plot(1:1:T,i_ea);
axis([0 40 -0.1 0.1])
title('Nominal interest rate ');
%% back out the price level and nominal money stock from monetary shock
T = 40;
lnP_em(1) = pi_em(1);
M_em(1) = m_em(1)+lnP_em(1);
for ii=2:T
    lnP_em(ii) = lnP_em(ii-1) + pi_em(ii); %in log, intial price level =1
    M_em(ii) = m_em(ii)+lnP_em(ii); %in log
end
figure;
subplot(2,2,1);
plot(1:1:T,lnP_em);
title('Price level ');
subplot(2,2,2);
plot(1:1:T,M_em);
axis([0 40 -0.1 0.1])
title('Nominal Money');
subplot(2,2,3);
plot(1:1:T,i_em);
axis([0 40 -0.1 0.1])
title('Nominal interest rate ');
subplot(2,2,4);
plot(1:1:T,yf_em);
axis([0 40 -0.1 0.1])
title('flexible output ');