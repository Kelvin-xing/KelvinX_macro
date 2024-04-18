clear all;
clc;

global rho alpha eta beta delta W r smin smax fspace cr egrid P theta lbar Pw

% fixed parameters
beta      = 0.92;            
alpha     = 2/3;
eta       = 0.85;               % decreasing returns
delta     = 0.06; 
Lbar      = 1;
lambdae  = 0.785;              % probability worker stays in e
lambdau  = 0.5;                % 1
spliorder = 1; 

% Targets

sdYdata   = 0.59; 
sYdata    = 1.31;
ac1Ydata  = 0.90;
ac3Ydata  = 0.87;
ac5Ydata  = 0.85;

DtoY      = 1.2;
EtoY      = 0.3;

%calibrated parameters 

rho    = 0.3;               % persistence transitory productivity
se     = 0.72;

W      = 0.812;
r      = 0.047;

theta1 = 0.15;
theta2 = 0.58; 
theta3 = 0.85;

% theta1 = 0; 
% theta2 = 0;
% theta3 = 0;
% W      = 0.752;

solve_workers;

disp('low theta')
theta = theta1; 
solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y1 = Ya;
K1 = Ka;
L1 = La; 
A1 = Aa;
C1 = n'*C;
n1 = n;
Debt1 = (K - state(:,1));

LK1 = L./K;
YLK1 = Y./(W*L + r*K); 
DK1  = (K - state(:,1))./K;

TFP1 = TFP;
TFP1best = TFPbest;
%save all relevant simulated micro-data for moments
simulate1
Ys1=Y;
Ls1=L;
Ks1=K;
es1=e;
rs1=rr;

disp('medium theta')
theta  = theta2; 
solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y2 = Ya;
K2 = Ka;
L2 = La; 
A2 = Aa;
C2 = n'*C;
n2 = n;
Debt2 = (K - state(:,1));

LK2 = L./K;
YLK2 = Y./(W*L + r*K); 
DK2  = (K - state(:,1))./K;

TFP2 = TFP;
TFP2best = TFPbest;
simulate1;
%save all relevant simulated micro-data for moments
Ys2=Y;
Ls2=L;
Ks2=K;
es2=e;
rs2=rr;

disp('high theta')
theta  = theta3; 
solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y3 = Ya;
K3 = Ka;
L3 = La; 
A3 = Aa;
C3 = n'*C;
n3 = n;
Debt3 = (K - state(:,1));

LK3 = L./K;
YLK3 = Y./(W*L + r*K); 
DK3  = (K - state(:,1))./K;

TFP3 = TFP;
TFP3best = TFPbest;
%save all relevant simulated micro-data for moments
simulate1;
Ys3=Y;
Ls3=L;
Ks3=K;
es3=e;
rs3=rr;

Ya = 1/3*(Y1 + Y2 + Y3); 
Ka = 1/3*(K1 + K2 + K3); 
La = 1/3*(L1 + L2 + L3); 
DtoYm = 1/3*(n1'*Debt1 + n2'*Debt2 + n3'*Debt3)/Ya;

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Ka,  Aa + Aw]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [La,  Lbar]);
fprintf('Debt to GDP                   = %9.3f  %9.3f \n',   [DtoYm,  DtoY]);

fprintf('\n');

D    = Ka - (Aa + Aw);                          % Net Foreign Debt Position
Ca   = Ya - delta*Ka - r*D; 
Ca2  = 1/3*(C1 + C2 + C3) + W*Lbar + r*Aw;      % add consumption of all agents

fprintf('\n');
disp('Model Implications')

fprintf('Output                        = %9.3f \n',          Ya/La);   % per worker in case L doesnt clear
fprintf('Consumption                   = %9.3f \n',          Ca/La);
fprintf('Investment                    = %9.3f \n',          delta*Ka/La);
fprintf('TFP  1                        = %9.3f \n',          TFP1);
fprintf('Misallocation Loss   1        = %9.3f \n',          log(TFP1best/TFP1));
fprintf('TFP  2                        = %9.3f \n',          TFP2);
fprintf('Misallocation Loss   2        = %9.3f \n',          log(TFP2best/TFP2));
fprintf('TFP  3                        = %9.3f \n',          TFP3);
fprintf('Misallocation Loss   3        = %9.3f \n',          log(TFP3best/TFP3));

state = gridmake(agrid, (1:1:k)');
TFP = Ya./(La.^alpha*Ka^(1-alpha))^eta;
TFPbest = (1/3*n1'*exp(egrid(state(:,2))) + 1/3*n2'*exp(egrid(state(:,2))) + 1/3*n3'*exp(egrid(state(:,2)))).^(1-eta);

fprintf('TFP  Total                    = %9.3f \n',          TFP);
fprintf('Misallocation Loss   Total    = %9.3f \n',          log(TFPbest/TFP));

weight = [n1; n2; n3];
weight = weight/sum(weight); 
DK = [DK1; DK2; DK3];

data = [weight, DK];
data = sortrows(data, 2); 
data(:,1) = cumsum(data(:,1));

[junk, i25] = min(abs(data(:,1) - 0.25)); 
[junk, i75] = min(abs(data(:,1) - 0.75)); 

fprintf('Ratio of D/K 75th to 25h per           = %9.3f \n',         data(i75,2)/data(i25,2));

sim_moments;

% 
% 
% 
% data = [weight, LK];
% data = sortrows(data, 2); 
% data(:,1) = cumsum(data(:,1));
% 
% [junk, i25] = min(abs(data(:,1) - 0.25)); 
% [junk, i75] = min(abs(data(:,1) - 0.75)); 
% 
% fprintf('\n');
% 
% fprintf('Ratio of LK 75th to 25h per            = %9.3f \n',         data(i75,2)/data(i25,2));
% 
% 
% YLK = [YLK1; YLK2; YLK3];
% 
% data = [weight, YLK];
% data = sortrows(data, 2); 
% data(:,1) = cumsum(data(:,1));
% 
% [junk, i25] = min(abs(data(:,1) - 0.25)); 
% [junk, i75] = min(abs(data(:,1) - 0.75)); 
% 
% fprintf('Ratio of Y/(wL+rK) 75th to 25h per    = %9.3f \n',         data(i75,2)/data(i25,2));
% 
