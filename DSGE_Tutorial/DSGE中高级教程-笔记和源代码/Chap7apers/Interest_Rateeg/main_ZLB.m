%% main file for the ZLB and interest rate peg
%please use dynare v4.4.3, version like 4.2.0 will not work;
clear all
dynare ZLB noclearall nolog
save oo_zlb_baseline ;

clear all;
dynare ZLB_four noclearall nolog;
save oo_zlb_four ;

clear all;
dynare ZLB_eight noclearall nolog;
save oo_zlb_eight ;

%% pluck out data 
clear all
%1st dimension: number of IRF simulated;
%2nd dimension:number of cases: H=0,4,8, totally 3 cases;
%3rd dimension: number of shocks: technology,government spending
%net export shock, totally 3 shocks; 
output   = zeros(40,3,3);
inflation = zeros(40,3,3);
int         = zeros(40,3,3);
realint   = zeros(40,3,3);
cons   = zeros(40,3,3);
labor= zeros(40,3,3);
wage = zeros(40,3,3);

load oo_zlb_baseline;
output(:,1,1) = y_ea;
output(:,1,2) =y_eg;
output(:,1,3) =y_ex;
inflation(:,1,1) = pi_ea;
inflation(:,1,2) = pi_eg;
inflation(:,1,3) = pi_ex;
int(:,1,1) = i_ea ;
int(:,1,2) =i_eg;
int(:,1,3) =i_ex;
realint(:,1,1)= r_ea;
realint(:,1,2) =r_eg;
realint(:,1,3) =r_ex;
cons(:,1,1)=c_ea;
cons(:,1,2)=c_eg;
cons(:,1,3)=c_ex;

labor(:,1,1)=n_ea;
labor(:,1,2)=n_eg;
labor(:,1,3)=n_ex;
wage(:,1,2)=w_eg;
wage(:,1,3)=w_ex;

load oo_zlb_four;
output(:,2,1) = y_ea;
output(:,2,2) =y_eg;
output(:,2,3) =y_ex;
inflation(:,2,1) = pi_ea;
inflation(:,2,2) = pi_eg;
inflation(:,2,3) = pi_ex;
int(:,2,1) = i_ea ;
int(:,2,2) =i_eg;
int(:,2,3) =i_ex;
realint(:,2,1)= r_ea;
realint(:,2,2) =r_eg;
realint(:,2,3) =r_ex;
cons(:,2,1)=c_ea;
cons(:,2,2)=c_eg;
cons(:,2,3)=c_ex;
labor(:,2,1)=n_ea;
labor(:,2,2)=n_eg;
labor(:,2,3)=n_ex;
wage(:,2,2)=w_eg;
wage(:,2,3)=w_ex;

load oo_zlb_eight;
output(:,3,1) = y_ea;
output(:,3,2) =y_eg;
output(:,3,3) =y_ex;
inflation(:,3,1) = pi_ea;
inflation(:,3,2) = pi_eg;
inflation(:,3,3) = pi_ex;
int(:,3,1) = i_ea ;
int(:,3,2) =i_eg;
int(:,3,3) =i_ex;
realint(:,3,1)= r_ea;
realint(:,3,2) =r_eg;
realint(:,3,3) =r_ex;
cons(:,3,1)=c_ea;
cons(:,3,2)=c_eg;
cons(:,3,3)=c_ex;
labor(:,3,1)=n_ea;
labor(:,3,2)=n_eg;
labor(:,3,3)=n_ex;
wage(:,3,2)=w_eg;
wage(:,3,3)=w_ex;

%% technology shock
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|-.|:')
T= 1:1:40;
subplot(2,3,1);
plot(T,output(:,:,1),'LineWidth',2);
title('output ');

subplot(2,3,2);
plot(T,inflation(:,:,1),'LineWidth',2);
title('inflation ');

subplot(2,3,3);
plot(T,int(:,:,1),'LineWidth',2);
title('int. rate.');

subplot(2,3,4);
plot(T,realint(:,:,1),'LineWidth',2);
title('real int. rate .');

subplot(2,3,5);
plot(T,cons(:,:,1),'LineWidth',2);
title('consumption');

subplot(2,3,6);
plot(T,labor(:,:,1),'LineWidth',2);
title('labor');

legend('H=0','H=4','H=8');
%% government shock
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.');
figure;
subplot(2,3,1);
plot(1:1:40,output(:,:,2),'LineWidth',2);
title('output');

subplot(2,3,2);
plot(1:1:40,cons(:,:,2),'LineWidth',2);
title('consumption');

subplot(2,3,3);
plot(1:1:40,int(:,:,2),'LineWidth',2);
title('int. rate');

subplot(2,3,4);
plot(1:1:40,realint(:,:,2),'LineWidth',2);
title('real int. rate');

subplot(2,3,5);
plot(1:1:40,inflation(:,:,2),'LineWidth',2);
title('inflation');

subplot(2,3,6);
plot(1:1:40,labor(:,:,2),'LineWidth',2);
title('labor');

legend('H=0','H=4','H=8');

%% net export shock
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.');
figure;
subplot(3,2,1);
plot(1:1:40,output(:,:,3),'LineWidth',2);
title('output');

subplot(3,2,2);
plot(1:1:40,cons(:,:,3),'LineWidth',2);
title('consumption');

subplot(3,2,3);
plot(1:1:40,int(:,:,3),'LineWidth',2);
title('nominal interest rate');

subplot(3,2,4);
plot(1:1:40,realint(:,:,3),'LineWidth',2);
title('real interest rate');
legend('H=0','H=4','H=8');

subplot(3,2,5)
plot(1:1:40,inflation(:,:,3),'LineWidth',2);
title('inflation');

subplot(3,2,6)
plot(1:1:40,labor(:,:,3),'LineWidth',2);
title('labor');
