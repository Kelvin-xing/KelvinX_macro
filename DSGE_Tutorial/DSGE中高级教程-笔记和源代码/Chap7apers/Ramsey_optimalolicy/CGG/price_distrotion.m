%we show how the price distortion evolution
clear all;
H=30;
pstar = zeros(2,H);
gstar = zeros(2,H);


epsilon=11;
theta=0.75;
pstar(1,1)=0.8;
pstar(2,1)=1.2;
gstar(1,1)=0.8;
gstar(2,1)=1.2;


for ii=2:H
    pstar(1,ii)=(1-theta+theta*pstar(1,ii-1)^(epsilon-1))^(1/((epsilon-1)));
    pstar(2,ii)=(1-theta+theta*pstar(2,ii-1)^(epsilon-1))^(1/((epsilon-1)));
    gstar(1,ii)=(1-theta+theta*gstar(1,ii-1)^(1-epsilon))^(1/((1-epsilon)));
    gstar(2,ii)=(1-theta+theta*gstar(2,ii-1)^(1-epsilon))^(1/((1-epsilon)));
end
T=1:1:H;
close all;
figure;
plot(T,pstar(1,:),'k*-',T,gstar(1,:),'bo-');
title('Initial value 0.8 for both pstar and gstar with \epsilon=11, \theta=0.75');
legend('pstar' ,'gstar')
figure;
plot(T,pstar(2,:),'k*-',T,gstar(2,:),'bo-');
legend('pstar' ,'gstar')
title('Initial value 1.2 for both pstar and gstar with \epsilon=11,\theta=0.75');