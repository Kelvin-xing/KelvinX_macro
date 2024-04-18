clear; 
clc;

eta   = 0.85; 
alpha = 2/3;
delta = 0.06;
r     = 0.04; 

gamma = 1.25; 
share = 0.5;

omega = [0.089
0.096
0.096
0.093
0.089
0.08
0.066
0.054
0.045
0.037
0.03
0.027
0.024
0.022
0.021
0.019
0.016
0.014
0.012
0.01
0.008
0.007
0.006
0.005
0.004
0.003
0.002
0.002
0.002
0.02];   % fraction in each age group (last includes 30+)

omega = omega/sum(omega); 

gY    = [0.511
0.162
0.159
0.169
0.106
0.127
0.061
0.059
0.101
0.113
0.127
0.067
0.040
0.096
0.042
0.022
0.039
0.065
0.096
-0.043
0.099
0.017
0.029
0.071
0.028
0.061
0.085
-0.076
0.078
0.059];


logY = cumsum([0; gY(1:end-1)]);

YK = - gamma*log((1:1:30)');

T1 = sum(omega.*exp(YK).^(-(1-alpha)*eta/(1-eta)))^(1-alpha*eta);
T2 = sum(omega.*exp(YK).^((alpha*eta-1)/(1-eta)))^((1-alpha)*eta);

TFPloss = -log(T1/T2);

disp('TFPloss')
disp(TFPloss)

disp('MPK 1-5 vs 11+')

MPK = exp(YK);

disp(log(omega(1:5)'*MPK(1:5)/sum(omega(1:5))/(omega(11:end)'*MPK(11:end)/sum(omega(11:end)))))



K   = MPK;
Y   = MPK; 
APK = MPK;

K(1) = 10; 
F = K(1)*share/(1-share);
Y(1) = MPK(1)*K(1);

Y    = exp(logY).*Y(1);
K    = Y.*MPK.^(-1);

APK = Y./(K+F);

disp('APK 1-5 vs 11+')
disp(log(omega(1:5)'*APK(1:5)/sum(omega(1:5))/(omega(11:end)'*APK(11:end)/sum(omega(11:end)))))


% Back out r

R = MPK./MPK(30)*(r+delta) - delta;

disp('shadow r  1-5 ')
disp(omega(1:5)'*R(1:5)/sum(omega(1:5)));

disp('shadow r  11+ ')
disp(omega(11:end)'*R(11:end)/sum(omega(11:end)));