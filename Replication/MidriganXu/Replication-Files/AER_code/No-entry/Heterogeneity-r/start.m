clear all;
clc;

global rho alpha eta beta delta W r smin smax fspace cr egrid P theta lbar Pw

beta      = 0.933;            

alpha     = 2/3;
eta       = 0.85;               % decreasing returns
delta     = 0.06; 
Lbar      = 1;

spliorder = 1; 

% Targets

sdYdata   = 0.89; 
sYdata    = 1.45;
ac1Ydata  = 0.80;
ac3Ydata  = 0.70;
ac5Ydata  = 0.65;

DtoY      = 0.7;

rho    = 0.50;               % persistence transitory productivity
se     = 0.90;

W      = 0.86;

rL     = 0.05;
rH1    = 0.15;
rH2    = (rH1+rL)/2;

r = rL; 
r1 = rL;                       % first group, guys that borrow at high interest rates
r2 = rL; 

solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y1 = Ya;
K1 = Ka;
L1 = La; 
A1 = Aa;
C1 = n'*C;
n1 = n;
Debt1 = (K - state(:,1));

LK1  = L./K;
YLK1 = Y./(W*L + r*K); 
DK1  = (K - state(:,1))./K;

TFP1 = TFP;
TFP1best = TFPbest;

simulate1

Ys1=Y;
Ls1=L;
Ks1=K;
es1=e;
rs1=rr;

%%%%%%%%%%%%%


r1 = rL;     % second group, both borrow at low interest rate
r2 = rH2; 


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


simulate1

Ys2=Y;
Ls2=L;
Ks2=K;
es2=e;
rs2=rr;

%%%%%%%%%%%%%


r1 = rL;     % second group, both borrow at low interest rate
r2 = rH1; 


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


simulate1

Ys3=Y;
Ls3=L;
Ks3=K;
es3=e;
rs3=rr;

%%%%%%%%%%%%%%%%%


Ya = 1/3*(Y1 + Y2 + Y3); 
Ka = 1/3*(K1 + K2 + K3); 
La = 1/3*(L1 + L2 + L3); 
DtoYm = 1/3*(n1'*Debt1 + n2'*Debt2 + n3'*Debt3)/Ya;

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [La,     Lbar]);
fprintf('Debt to GDP                   = %9.3f  %9.3f \n',   [DtoYm,  DtoY]);

fprintf('\n');

fprintf('\n');
disp('Model Implications')


fprintf('TFP  1                        = %9.3f \n',          TFP1);
fprintf('Misallocation Loss   1        = %9.3f \n',          log(TFP1best/TFP1));
fprintf('TFP  2                        = %9.3f \n',          TFP2);
fprintf('Misallocation Loss   2        = %9.3f \n',          log(TFP2best/TFP2));
fprintf('TFP  2                        = %9.3f \n',          TFP3);
fprintf('Misallocation Loss   3        = %9.3f \n',          log(TFP3best/TFP3));

state = gridmake(agrid, (1:1:k)');

TFP     = Ya./(La.^alpha*Ka^(1-alpha))^eta;
TFPbest = (1/3*n1'*exp(egrid(state(:,2))) + 1/3*n2'*exp(egrid(state(:,2))) + 1/3*n3'*exp(egrid(state(:,2)))).^(1-eta);

fprintf('TFP  Total                    = %9.3f \n',          TFP);
fprintf('Misallocation Loss   Total    = %9.3f \n',          log(TFPbest/TFP));

sim_moments
