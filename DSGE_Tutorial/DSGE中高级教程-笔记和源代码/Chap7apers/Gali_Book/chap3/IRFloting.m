%this m script should be run after runing chap3.mod
%% IRF plotting: monetary shock
%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.1, p. 53
%----------------------------------------------------------------
close all
T=1:1:15;
figure(1)
subplot(3,2,1)
plot(T,y_gap_eps_nu,'-o');
title('Output Gap');

subplot(3,2,2)
plot(T,pi_ann_eps_nu,'-o');
title('Inflation');

subplot(3,2,3)
plot(T,R_ann_eps_nu,'-o');
title('Nominal Rate');

subplot(3,2,4)
plot(T,r_ann_eps_nu,'-o');
title('Real Rate');

subplot(3,2,5)
plot(T,m_growth_ann_eps_nu,'-o');
title('Money Growth');

subplot(3,2,6)
plot(T,nu_eps_nu,'-o');
title('\nu');

%% IRF plotting: technology shock
%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.2, p. 55
%----------------------------------------------------------------
close all
T=1:1:15;
figure(1)
subplot(4,2,1)
plot(T,y_gap_eps_a,'-o');
title('Output Gap');

subplot(4,2,2)
plot(T,pi_ann_eps_a,'-o');
title('Inflation');

subplot(4,2,3)
plot(T,y_eps_a,'-o');
title('Output');

subplot(4,2,4)
plot(T,n_eps_a,'-o');
title('Empolyment');

subplot(4,2,5)
plot(T,R_ann_eps_a,'-o');
title('Nominal Rate');

subplot(4,2,6)
plot(T,r_ann_eps_a,'-o');
title('Real Rate');

subplot(4,2,7)
plot(T,m_growth_ann_eps_a,'-o');
title('Money Growth');

subplot(4,2,8)
plot(T,a_eps_a,'-o');
title('Technology Shock');
